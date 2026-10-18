"""Tests for the real-data adapter (poa.services.realdata_import)."""
import json

from poa.core.parsers import bepipred
from poa.services import realdata_import


def _write_bepipred2_json(path, seq, pred):
    """Minimal BepiPred-2.0-style JSON: single generic 'Sequence' antigen, one score per residue."""
    data = {
        "info": {"failedjobs": 0, "size": 1},
        "antigens": {"Sequence": {"AA": list(seq), "PRED": pred}},
    }
    path.write_text(json.dumps(data))
    return str(path)


def test_specie_from_filename():
    assert realdata_import.specie_from_filename("/x/bepipred_denv1.json") == "DENV1"
    assert realdata_import.specie_from_filename("bepipred_chikv.json") == "CHIKV"


def test_prepare_rebuilds_reference_and_conforming_key(tmp_path):
    seq = "MRCVGIGN"
    pred = [0.1, 0.9, 0.9, 0.9, 0.1, 0.1, 0.9, 0.9]
    src = _write_bepipred2_json(tmp_path / "bepipred_denv1.json", seq, pred)

    prep = realdata_import.prepare_bepipred2(src, "DENV1", str(tmp_path / "out"), protein="E")

    assert prep.header == "E_DENV1_ref"
    assert prep.sequence == seq
    assert prep.has_x is False

    # reference FASTA is the exact embedded sequence
    ref = (tmp_path / "out" / "DENV1_ref.fasta").read_text().splitlines()
    assert ref[0] == ">E_DENV1_ref"
    assert ref[1] == seq

    # the rewritten JSON now drives correct species/protein extraction downstream
    b2 = bepipred.bp2_JsonAnalysis(prep.b2_json)
    df = bepipred.bp2_AntigenEpitopes(b2)
    assert set(df["Specie"]) == {"DENV1"}
    assert set(df["Protein"]) == {"E"}
    # PRED>0.5 residues form epitopes: positions 2-4 (RCV..) and 7-8
    assert not df.empty


def test_has_x_flag_is_set(tmp_path):
    src = _write_bepipred2_json(tmp_path / "bepipred_denv3.json", "MRCXGIG", [0.9] * 7)
    prep = realdata_import.prepare_bepipred2(src, "DENV3", str(tmp_path / "out"))
    assert prep.has_x is True


def test_best_world_match_picks_highest_identity(tmp_path):
    world = tmp_path / "world"
    world.mkdir()
    (world / "alpha.fasta").write_text(">a\nMRCVGIGNAAAA\n")   # identical prefix
    (world / "beta.fasta").write_text(">b\nWWWWWWWWWWWW\n")     # no match
    path, ident = realdata_import.best_world_match("MRCVGIGN", str(world))
    assert path.endswith("alpha.fasta")
    assert ident == 100.0


def test_best_world_match_tiebreak_prefers_closest_length(tmp_path):
    """Two near-identical hits: prefer the one whose length matches the reference (single protein)."""
    ref = "MRCVGIGNAA"  # len 10
    world = tmp_path / "world"
    world.mkdir()
    # single-protein set: same length, 100% identity
    (world / "single.fasta").write_text(">s\nMRCVGIGNAA\n")
    # polyprotein containing the reference: also 100% identity but much longer
    (world / "poly.fasta").write_text(">p\n" + "MRCVGIGNAA" + "K" * 40 + "\n")
    path, ident = realdata_import.best_world_match(ref, str(world))
    assert path.endswith("single.fasta")
    assert ident == 100.0
