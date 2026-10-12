"""Tests for POA0: the database search and the files it writes.

All network calls are stubbed. The payload shapes below are trimmed copies of real responses
(UniProt A0A0D4BSW3, PDBe 6h3x), including the parts that are easy to get wrong.
"""
import json

import pytest

from poa.core import dbsearch
from poa.services import database_client as db

# PDBe reports two numberings for the same residue. For 6H3X the author number is 482 where the
# entity index is 1, and the author number of the *last* residue is null. Reading the entity
# index as if it were the file's numbering shifts every epitope by 481 residues.
SIFTS_6H3X = {
    "6h3x": {"UniProt": {"A0A0D4BSW3": {"identifier": "A0A0D4BSW3_9VIRU", "mappings": [{
        "entity_id": 1, "chain_id": "A", "struct_asym_id": "A",
        "unp_start": 482, "unp_end": 703,
        "start": {"author_residue_number": 482, "author_insertion_code": "", "residue_number": 1},
        "end": {"author_residue_number": None, "author_insertion_code": "", "residue_number": 222},
        "identity": 1.0, "coverage": 0.156,
    }]}}}
}

BEST_STRUCTURES = {"Q32ZE1": [
    {"pdb_id": "5u0b", "chain_id": "A", "unp_start": 2517, "unp_end": 3419,
     "coverage": 0.264, "resolution": 3.0},
    {"pdb_id": "5u0b", "chain_id": "B", "unp_start": 2517, "unp_end": 3419,
     "coverage": 0.264, "resolution": 3.0},
    {"pdb_id": "5ire", "chain_id": "A", "unp_start": 291, "unp_end": 790,
     "coverage": 0.146, "resolution": 3.8},
]}

PDB_HEADER_4C23 = """\
HEADER    TRANSFERASE                             16-AUG-13   4C23
TITLE     L-FUCULOSE KINASE
COMPND   2 MOLECULE: L-FUCULOSE KINASE FUCK;
COMPND   3 CHAIN: A, B;
SOURCE   2 ORGANISM_SCIENTIFIC: STREPTOCOCCUS PNEUMONIAE;
EXPDTA    X-RAY DIFFRACTION
REMARK   2 RESOLUTION.    2.00 ANGSTROMS.
SEQRES   1 A  473  LEU VAL PRO ARG GLY SER HIS GLU ALA VAL LEU ALA ILE
SEQRES   1 B  473  LEU VAL PRO ARG GLY SER HIS GLU ALA VAL LEU ALA ILE
ATOM      1  CA  LEU A  -5      0.000   0.000   0.000  1.00  0.00           C
"""


# --------------------------------------------------------------------------- header tokens
@pytest.mark.parametrize("raw, expected", [
    ("Envelope protein E", "EnvelopeproteinE"),
    ("RNA-directed RNA polymerase", "RNAdirectedRNApolymerase"),
    ("Gc_head", "Gchead"),          # the underscore would create a phantom field
    ("", "SEQ"),
])
def test_sanitize_token_removes_what_splits_a_header(raw, expected):
    """
    Headers and epitope names are split on ``_``, so a token containing one shifts every field
    after it: ``poa.core.topology.sameVirus`` would read the species out of the protein's name.
    """
    assert dbsearch.sanitize_token(raw) == expected


def test_species_token_from_organism_name():
    assert dbsearch.species_token("Oropouche virus") == "OROPOUCHE"
    assert dbsearch.species_token("") == "ORG"


def test_written_header_has_the_three_fields_the_pipeline_reads(tmp_path):
    from poa.core.topology import sameVirus

    path = dbsearch.write_fasta(str(tmp_path / "x.fasta"), "GcHead_OROV_A0A0D4BSW3", "MKWV")
    header = open(path, encoding="utf-8").readline().strip().lstrip(">")

    fields = header.split("_")
    assert len(fields) == 3 and fields[0] == "GcHead" and fields[1] == "OROV"
    # This is the comparison POA2 actually performs when matching epitopes to a reference.
    assert sameVirus("OROV", header) is True
    assert sameVirus("DENV1", header) is False


def test_write_fasta_wraps_and_round_trips(tmp_path):
    from Bio import SeqIO

    sequence = "ACDEFGHIKLMNPQRSTVWY" * 8        # 160 aa, longer than one line
    path = dbsearch.write_fasta(str(tmp_path / "s.fasta"), "P_S_ACC", sequence)

    record = next(SeqIO.parse(path, "fasta"))
    assert str(record.seq) == sequence and record.id == "P_S_ACC"


# --------------------------------------------------------------------------- SIFTS parsing
def test_sifts_mapping_uses_author_numbering_not_the_entity_index(monkeypatch):
    monkeypatch.setattr(db, "_get_json", lambda *a, **k: SIFTS_6H3X)

    mapping = db.sifts_mappings("6H3X", accession="A0A0D4BSW3")[0]

    assert (mapping.unp_start, mapping.unp_end) == (482, 703)
    assert mapping.pdb_start == 482            # author numbering, what the file says
    assert mapping.entity_start == 1           # entity index, what it does NOT say
    assert mapping.offset == 0                 # taking the entity index would give -481


def test_sifts_derives_the_end_when_pdbe_omits_it(monkeypatch):
    monkeypatch.setattr(db, "_get_json", lambda *a, **k: SIFTS_6H3X)

    mapping = db.sifts_mappings("6H3X")[0]

    assert mapping.derived_end is True
    assert mapping.pdb_end == 703              # 482 + (703 - 482)


def test_sifts_filters_by_accession(monkeypatch):
    monkeypatch.setattr(db, "_get_json", lambda *a, **k: SIFTS_6H3X)

    assert db.sifts_mappings("6H3X", accession="P12345") == []


def test_sifts_absent_is_empty_not_an_error(monkeypatch):
    def unavailable(*a, **k):
        raise db.ServiceUnavailable("404")

    monkeypatch.setattr(db, "_get_json", unavailable)
    assert db.sifts_mappings("1ABC") == []


# --------------------------------------------------------------------------- listings
def test_list_structures_groups_chains_of_one_entry(monkeypatch):
    """
    One request for the whole accession, not one per structure.

    Walking the entry's PDB cross-references costs a request each; Zika's polyprotein has 281 of
    them, which turned a listing into minutes of sequential HTTP.
    """
    monkeypatch.setattr(db, "_get_json", lambda *a, **k: BEST_STRUCTURES)

    structures = db.list_structures("Q32ZE1")

    assert [s.pdb_id for s in structures] == ["5U0B", "5IRE"]      # best coverage first
    assert structures[0].chains == ["A", "B"]                      # one entry, two chains
    assert (structures[0].unp_start, structures[0].unp_end) == (2517, 3419)


def test_search_proteins_ranks_entries_with_a_structure_first(monkeypatch):
    payload = {"results": [
        {"primaryAccession": "A1", "entryType": "UniProtKB unreviewed (TrEMBL)",
         "sequence": {"length": 2252},
         "proteinDescription": {"submissionNames": [{"fullName": {"value": "RdRp"}}]}},
        {"primaryAccession": "A2", "entryType": "UniProtKB unreviewed (TrEMBL)",
         "sequence": {"length": 1420},
         "proteinDescription": {"submissionNames": [{"fullName": {"value": "Polyprotein"}}]},
         "uniProtKBCrossReferences": [{"database": "PDB", "id": "6H3X"}]},
    ]}
    monkeypatch.setattr(db, "_get_json", lambda *a, **k: payload)

    proteins = db.search_proteins(118655, limit=5)

    # The shorter entry wins because it has a structure: that is what makes a target usable.
    assert proteins[0].accession == "A2" and proteins[0].pdb_ids == ["6H3X"]
    assert proteins[0].database == "TrEMBL"


def test_search_proteins_overfetches_so_ranking_can_see_past_the_limit(monkeypatch):
    seen = {}

    def capture(url, params=None, **k):
        seen.update(params or {})
        return {"results": []}

    monkeypatch.setattr(db, "_get_json", capture)
    db.search_proteins(118655, limit=4)

    # Asking for exactly `limit` rows would let the one entry with a structure fall off the page
    # before the sort ever sees it.
    assert seen["size"] > 4
    assert seen["query"] == "organism_id:118655"


def test_include_descendants_switches_the_query_field(monkeypatch):
    seen = {}

    def capture(url, params=None, **k):
        seen.update(params or {})
        return {"results": []}

    monkeypatch.setattr(db, "_get_json", capture)

    db.search_proteins(11320, include_descendants=True, limit=3)
    # organism_id matches one taxon; taxonomy_id includes strains and isolates, which for a virus
    # can be an order of magnitude more entries.
    assert seen["query"] == "taxonomy_id:11320"


def test_regions_merge_sources_instead_of_dropping_one(monkeypatch):
    """A TrEMBL Domain is derived from Pfam, so both report the same span; neither is discarded."""
    entry = {"sequence": {"length": 1420}, "features": [
        {"type": "Domain", "description": "Bunyavirus glycoprotein G1",
         "location": {"start": {"value": 514}, "end": {"value": 1359}}}]}
    monkeypatch.setattr(db, "_pfam_regions", lambda acc: [
        db.Region(name="Bunyavirus glycoprotein G1", start=514, end=1359, source="Pfam")])
    monkeypatch.setattr(db, "list_structures", lambda acc, limit=10: [])

    regions = db.list_regions("A0A0D4BSW3", entry=entry)

    assert len(regions) == 1
    assert set(regions[0].source.split("+")) == {"Domain", "Pfam"}


def test_whole_length_chain_is_not_offered_as_a_region(monkeypatch):
    """A SwissProt polyprotein annotates itself as a Chain; that is the protein, not a part."""
    entry = {"sequence": {"length": 3419}, "features": [
        {"type": "Chain", "description": "Genome polyprotein",
         "location": {"start": {"value": 1}, "end": {"value": 3419}}},
        {"type": "Chain", "description": "Envelope protein E",
         "location": {"start": {"value": 291}, "end": {"value": 790}}}]}
    monkeypatch.setattr(db, "_pfam_regions", lambda acc: [])
    monkeypatch.setattr(db, "list_structures", lambda acc, limit=10: [])

    regions = db.list_regions("Q32ZE1", entry=entry)

    assert [r.name for r in regions] == ["Envelope protein E"]
    assert regions[0].length == 500


# --------------------------------------------------------------------------- structure header
def test_structure_summary_reads_what_the_file_actually_is():
    """
    A PDB id carries no meaning on sight, so a file can stand in for an unrelated protein with
    no visible symptom. 4C23 was used as the structure for a dengue analysis.
    """
    summary = db.structure_summary(PDB_HEADER_4C23)

    assert summary["molecule"] == "L-FUCULOSE KINASE FUCK"
    assert summary["organism"] == "STREPTOCOCCUS PNEUMONIAE"
    assert summary["method"] == "X-RAY DIFFRACTION"
    assert summary["chains"] == "A,B"


# --------------------------------------------------------------------------- region choice
def test_resolve_region_accepts_an_explicit_span():
    region = dbsearch.resolve_region("A0A0D4BSW3", "482-703", entry={})

    assert (region.start, region.end, region.length) == (482, 703, 222)
    assert region.source == "manual"


def test_resolve_region_matches_a_name(monkeypatch):
    monkeypatch.setattr(db, "list_regions", lambda acc, entry=None: [
        db.Region(name="Envelope protein E", start=291, end=790, source="Chain")])

    assert dbsearch.resolve_region("Q32ZE1", "envelope protein e", entry={}).start == 291
    assert dbsearch.resolve_region("Q32ZE1", "envelope", entry={}).start == 291
    assert dbsearch.resolve_region("Q32ZE1", "capsid", entry={}) is None


# --------------------------------------------------------------------------- end to end
def test_run_poa0_writes_the_fasta_and_the_report(tmp_path, monkeypatch):
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 100          # 2000 aa
    monkeypatch.setattr(db, "get_entry", lambda acc: {
        "sequence": {"length": len(sequence), "value": sequence},
        "organism": {"scientificName": "Oropouche virus"},
        "proteinDescription": {"submissionNames": [{"fullName": {"value": "Envelopment polyprotein"}}]},
        "features": []})

    class Args:
        s, a, d = "", "A0A0D4BSW3", str(tmp_path)
        r, nome, sigla, estrutura = "482-703", "GcHead", "OROV", ""
        descendentes, limite, l = False, 25, ""

    report = dbsearch.run_poa0(Args())

    written = open(report["arquivos"]["fasta"], encoding="utf-8").read().splitlines()
    assert written[0] == ">GcHead_OROV_A0A0D4BSW3"
    assert "".join(written[1:]) == sequence[481:703]
    assert len("".join(written[1:])) == 222

    saved = json.load(open(report["arquivos"]["relatorio"], encoding="utf-8"))
    assert saved["regiao"]["start"] == 482 and saved["sequencia"]["tamanho"] == 222
