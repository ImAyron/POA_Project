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
    assert report.sources == [(1, 5)]


@pytest.mark.parametrize("sifts", [False, True])
def test_every_mapped_range_can_be_traced_to_its_epitope(sifts):
    """
    ``ranges`` leaves out what did not map, so its positions do not line up with the input.

    A viewer that names the stretches it colors needs to know which epitope each one is. With an
    unmapped epitope in the middle, the second range is the *third* epitope, and ``sources`` is
    the only thing that says so.
    """
    structure = make_pdb(REF[:12], first=482)                 # covers 1-12 of the sequence
    wanted = [(1, 4), (15, 19), (8, 11)]                      # the middle one is not in the file
    extra = dict(sifts_offset=481, sifts_span=(1, 12)) if sifts else {}

    report = sm.map_epitopes(structure, REF, wanted, **extra)

    assert report.ok is True
    assert report.sources == [(1, 4), (8, 11)]
    assert dict(zip(report.sources, report.ranges)) == {(1, 4): (482, 485), (8, 11): (489, 492)}
    assert report.unmapped == [(15, 19)]


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


# --------------------------------------------------------------------------- which chain
def test_the_chain_that_is_the_protein_is_found_wherever_it_is():
    """
    In a complex the antigen is often not the first chain.

    With no chain named, the first one was compared, found to be another protein, and the whole
    structure refused — while the right chain sat further down the same file.
    """
    complex_ = make_pdb(OTHER, chain="H") + make_pdb(REF, chain="B", first=482)

    assert sm.map_epitopes(complex_, REF, [(1, 5)]).ok is False          # the old behaviour
    report = sm.map_epitopes_best_chain(complex_, REF, [(1, 5)])

    assert report.ok is True
    assert report.chain == "B"
    assert report.ranges == [(482, 486)]


def test_the_first_of_equal_chains_is_kept():
    """Two copies map equally well; the choice must be stable, and file order is the stable one."""
    dimer = make_pdb(REF, chain="A") + make_pdb(REF, chain="B")

    assert sm.map_epitopes_best_chain(dimer, REF, [(1, 5)]).chain == "A"


def test_a_full_chain_beats_a_fragment_of_the_same_protein():
    both = make_pdb(REF[:12], chain="A") + make_pdb(REF, chain="B")

    assert sm.map_epitopes_best_chain(both, REF, [(1, 5)]).chain == "B"


def test_no_matching_chain_is_still_a_refusal_and_says_how_many_were_tried():
    unrelated = make_pdb(OTHER, chain="A") + make_pdb(OTHER[::-1], chain="B")

    report = sm.map_epitopes_best_chain(unrelated, REF, [(1, 5)])

    assert report.ok is False and report.ranges == []
    assert "Nenhuma das 2 cadeias" in report.warnings[0]


def test_a_file_with_no_protein_chain_is_refused_not_raised():
    assert sm.map_epitopes_best_chain("END\n", REF, [(1, 5)]).ok is False


def test_residues_by_chain_agrees_with_reading_each_chain():
    text = make_pdb(REF, chain="A", first=-5) + make_pdb(OTHER, chain="B", skip=(3, 4))

    by_chain = sm.residues_by_chain(text)

    assert list(by_chain) == ["A", "B"]
    for chain in by_chain:
        assert by_chain[chain] == sm.chain_residues(text, chain)


def test_copies_are_chains_that_call_the_same_residues_by_the_same_numbers():
    """
    A dimer has the protein twice; a disordered loop in one copy does not make it another protein.

    What does disqualify a chain is numbering: the ranges found for one chain are residue
    numbers, and they only mean the same residues on a chain numbered the same way.
    """
    text = (make_pdb(REF, chain="A", first=482)
            + make_pdb(REF, chain="B", first=482, skip=(7, 8))      # same, one loop missing
            + make_pdb(REF, chain="C", first=1)                     # same protein, other numbers
            + make_pdb(OTHER, chain="D", first=482))                # same numbers, other protein

    assert sm.chain_copies(text, "A") == ["B"]
    assert sm.chain_copies(text, "B") == ["A"]
    assert sm.chain_copies(text, "Z") == []


# --------------------------------------------------------------------------- mmCIF
CIF_COLUMNS = ("group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
               "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
               "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
               "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num")


def make_cif(sequence, *, chain="A", label_chain="A", first=1, skip=(), models=1, body_only=False):
    """
    The mmCIF counterpart of :func:`make_pdb`: one CA row per residue.

    ``label_seq_id`` always counts from 1 while ``auth_seq_id`` starts at ``first`` — the two
    numberings a real file carries for the same residue.
    """
    rows, serial = [], 0
    for model in range(1, models + 1):
        for index, letter in enumerate(sequence, start=1):
            if index in skip:
                continue
            serial += 1
            name = THREE[letter]
            rows.append(f"ATOM {serial} C CA . {name} {label_chain} 1 {index} ? {index:.3f} 0.000 "
                        f"0.000 1.00 0.00 ? {first + index - 1} {name} {chain} CA {model}")
    body = "\n".join(rows) + "\n"
    if body_only:
        return body
    header = ("data_TEST\n#\n_struct.title 'A test structure'\n"
              "_exptl.method 'X-RAY DIFFRACTION'\n_refine.ls_d_res_high 2.10\n#\n"
              "loop_\n" + "\n".join(f"_atom_site.{column}" for column in CIF_COLUMNS) + "\n")
    return header + body + "#\n"


def test_the_format_is_read_from_the_content():
    assert sm.structure_format(make_pdb(REF)) == "pdb"
    assert sm.structure_format(make_cif(REF)) == "cif"
    assert sm.structure_format("") == "pdb"


def test_an_mmcif_file_reads_the_same_as_its_pdb_counterpart():
    """
    The numbering that counts is the authors', in both formats.

    An mmCIF row carries ``label_seq_id`` (1, 2, 3… within the entity) beside ``auth_seq_id``
    (what the PDB-format file of the same entry says). For 6H3X they are 1 and 482. Reading the
    wrong one would place every epitope 481 residues away — and the viewer, which reads the
    ``auth_`` fields, would color whatever it was told.
    """
    cif = make_cif(REF, chain="B", label_chain="A", first=482, skip=(7,))
    pdb = make_pdb(REF, chain="B", first=482, skip=(7,))

    assert sm.chain_ids(cif) == sm.chain_ids(pdb) == ["B"]
    assert sm.chain_residues(cif) == sm.chain_residues(pdb)
    assert sm.residues_by_chain(cif) == sm.residues_by_chain(pdb)
    assert sm.chain_residues(cif)[0].number == 482

    from_cif = sm.map_epitopes_best_chain(cif, REF, [(1, 5), (10, 15)])
    from_pdb = sm.map_epitopes_best_chain(pdb, REF, [(1, 5), (10, 15)])
    assert from_cif.ok and from_cif.ranges == from_pdb.ranges == [(482, 486), (491, 496)]
    assert from_cif.partial == from_pdb.partial


def test_only_the_first_model_of_an_mmcif_file_is_read():
    assert len(sm.chain_residues(make_cif(REF, models=3))) == len(REF)


def test_mmcif_chains_with_long_identifiers_are_kept_whole():
    """Entries too large for the PDB format are the ones with chains named 'AA', 'a1', …"""
    cif = make_cif(REF, chain="AA") .rstrip("#\n") + "\n" + make_cif(OTHER, chain="a1", body_only=True)

    assert sm.chain_ids(cif) == ["AA", "a1"]
    assert sm.map_epitopes_best_chain(cif, REF, [(1, 5)]).chain == "AA"


def test_an_epitope_with_no_residue_in_the_file_is_unmapped_under_sifts_too():
    """
    SIFTS says where a residue *would* be numbered; that is not the residue being there.

    With the span covering the whole reference, an epitope outside the stretch the crystal
    resolves came back as mapped and merely "partial" — a range of residue numbers with nothing
    in it, counted among the highlighted ones. The alignment path already called it unmapped.
    """
    structure = make_pdb(REF[:10], first=482)                 # resolves 1-10 of the sequence

    report = sm.map_epitopes(structure, REF, [(1, 5), (8, 12), (15, 19)],
                             sifts_offset=481, sifts_span=(1, len(REF)))

    assert report.ok is True
    assert report.sources == [(1, 5), (8, 12)]
    assert report.partial == [(8, 12)]                        # three of its five residues exist
    assert report.unmapped == [(15, 19)]
