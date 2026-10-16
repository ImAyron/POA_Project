"""Tests for mapping epitope positions onto a structure's residue numbering.

The case these guard against is the one that fails silently: a structure of a *different*
protein whose residue numbering happens to span the epitope positions. 3Dmol.js colors whatever
numbers it is handed, so the only symptom is a plausible, meaningless figure.
"""
import pytest

from poa.core import structure_map as sm

# A 20-residue sequence with no internal repeats, so an alignment has one obvious answer.
REF = "MKWVTFISLLFLFSSAYSRG"
# Something unrelated of the same length, to stand in for "a structure of another protein".
OTHER = "PQDGHCNEYIPQDGHCNEYI"

THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
         "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
         "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def make_pdb(sequence, *, chain="A", first=1, skip=()):
    """
    Build a minimal PDB with one CA record per residue.

    Args:
        sequence: one-letter sequence.
        chain: chain identifier.
        first: residue number of the first residue — the knob that reproduces a structure
            numbered from something other than 1.
        skip: 1-based offsets into ``sequence`` to leave out, standing in for residues that are
            present in the construct but disordered in the crystal.
    """
    lines = []
    for index, letter in enumerate(sequence, start=1):
        if index in skip:
            continue
        number = first + index - 1
        lines.append(
            f"ATOM  {index:>5}  CA  {THREE[letter]} {chain}{number:>4}    "
            f"{index:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
    lines.append("END")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- parsing
def test_chain_residues_reads_numbering_from_the_file():
    residues = sm.chain_residues(make_pdb(REF, first=482))
    assert len(residues) == len(REF)
    assert residues[0].number == 482 and residues[0].letter == "M"
    assert residues[-1].number == 482 + len(REF) - 1


def test_chain_ids_and_default_chain():
    pdb = make_pdb(REF, chain="B")
    assert sm.chain_ids(pdb) == ["B"]
    # With no chain requested the first one found is used, rather than assuming "A" exists.
    assert sm.chain_sequence(pdb) == REF


def test_negative_residue_numbers_are_read():
    """Expression tags push numbering below 1 — PDB 4C23's chain A starts at -5."""
    residues = sm.chain_residues(make_pdb(REF, first=-5))
    assert residues[0].number == -5


# --------------------------------------------------------------------------- the real guard
def test_structure_of_another_protein_is_refused():
    """
    The motivating failure: 19 DENV epitopes highlighted on L-fuculose kinase (PDB 4C23).

    None of the epitope sequences occur in that structure, yet every position falls inside its
    numbering, so the viewer painted all of them with no error.
    """
    report = sm.map_epitopes(make_pdb(OTHER), REF, [(1, 5), (10, 15)])

    assert report.ok is False
    assert report.ranges == []                  # nothing is offered for highlighting
    assert report.identity < sm.DEFAULT_MIN_IDENTITY
    assert "outra proteina" in " ".join(report.warnings)


def test_matching_structure_maps_and_keeps_the_files_numbering():
    report = sm.map_epitopes(make_pdb(REF, first=482), REF, [(1, 5), (10, 15)])

    assert report.ok is True
    assert report.identity == pytest.approx(1.0)
    # Positions are reported in the numbering the viewer will use, not the sequence's.
    assert report.ranges == [(482, 486), (491, 496)]


def test_epitope_partly_missing_from_the_crystal_is_mapped_but_flagged():
    """
    An epitope whose middle is disordered is still placed correctly — the viewer draws the
    residues that exist. What the picture cannot show is that it is drawing only half of it,
    so the report says so rather than passing it off as a full match.
    """
    report = sm.map_epitopes(make_pdb(REF, skip=(11, 12, 13)), REF, [(1, 5), (10, 15)])

    assert report.ok is True
    assert report.ranges == [(1, 5), (10, 15)]
    assert (10, 15) in report.partial and (1, 5) not in report.partial
    assert "parcialmente ausentes" in " ".join(report.warnings)


def test_epitope_with_a_missing_end_is_clipped_to_what_was_observed():
    """
    A missing end shortens the range instead of extending it by arithmetic.

    Without SIFTS there is no way to know the residue number of a residue the file does not
    contain: numbering can jump or carry insertion codes, so "last observed + 1" is a guess.
    Clipping uses only positions that genuinely aligned, and the epitope is flagged as partial.

    The removed residue is position 16 ("A"), the only alanine in ``REF``. Dropping one of a
    repeated pair instead would make this untestable: removing the second S of "SS" leaves an
    alignment that cannot tell which of the two went missing, and either answer is as good.
    """
    assert REF.count("A") == 1 and REF[15] == "A"
    report = sm.map_epitopes(make_pdb(REF, skip=(16,)), REF, [(12, 16)])

    assert report.ok is True
    assert report.ranges == [(12, 15)]          # clipped, not extended to 16
    assert (12, 16) in report.partial


def test_epitope_entirely_absent_from_the_structure_is_unmapped():
    """With nothing observed there is no anchor at all, so the epitope is reported as missing."""
    report = sm.map_epitopes(make_pdb(REF[:10]), REF, [(15, 20)])

    assert report.ok is True
    assert (15, 20) in report.unmapped
    assert report.ranges == []


def test_missing_chain_is_an_explicit_failure():
    report = sm.map_epitopes(make_pdb(REF, chain="A"), REF, [(1, 5)], chain="Z")

    assert report.ok is False
    assert "Z" in report.warnings[0] and "A" in report.warnings[0]


def test_identity_below_full_match_is_flagged_but_usable():
    """A structure of a close variant is still worth showing, with the difference stated."""
    variant = "MKWVTFISLLFLFSSAYSKG".replace("SKG", "SKA")    # 2 substitutions in 20
    report = sm.map_epitopes(make_pdb(variant), REF, [(1, 10)])

    assert report.ok is True
    assert 0.70 <= report.identity < 1.0
    assert any("variante" in w for w in report.warnings)


# --------------------------------------------------------------------------- SIFTS path
def test_sifts_offset_shifts_positions_without_aligning():
    report = sm.map_epitopes(make_pdb(REF, first=482), REF, [(1, 5)],
                             sifts_offset=481, sifts_span=(1, len(REF)))

    assert report.method == "sifts"
    assert report.ok is True
    assert report.ranges == [(482, 486)]


def test_wrong_sifts_offset_is_caught_by_checking_the_residues():
    """
    PDBe reports two numberings and they differ: for 6H3X, ``residue_number`` is 1 where
    ``author_residue_number`` is 482. Taking one for the other shifts every epitope by 481, so
    the offset is verified against the actual residues instead of being trusted.
    """
    report = sm.map_epitopes(make_pdb(REF, first=482), REF, [(1, 5)],
                             sifts_offset=0, sifts_span=(1, len(REF)))

    assert report.ok is False
    assert report.ranges == []


def test_epitope_outside_the_sifts_span_is_unmapped():
    report = sm.map_epitopes(make_pdb(REF, first=482), REF, [(1, 5), (18, 20)],
                             sifts_offset=481, sifts_span=(1, 15))

    assert report.ok is True
    assert (18, 20) in report.unmapped
    assert report.ranges == [(482, 486)]


# --------------------------------------------------------------------------- alignment
def test_align_sequences_maps_a_subsequence_to_its_place():
    """A structure covering the middle of a longer sequence must land where it belongs."""
    index_map, identity, columns = sm.align_sequences(REF, REF[5:15])

    assert identity == pytest.approx(1.0)
    assert columns == 10
    assert index_map[5] == 0 and index_map[14] == 9


def test_empty_inputs_do_not_raise():
    assert sm.align_sequences("", "ABC") == ({}, 0.0, 0)
    report = sm.map_epitopes("END\n", REF, [(1, 5)])
    assert report.ok is False
