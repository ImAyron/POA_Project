"""Tests for the multi-species run (two or more organisms in the same POA1/Conservancy pass).

Covers the three failures that made a two-species run collapse into one:
  * the BepiPred-3.0 parser read the header as ``Specie_Protein_ID`` (swapped);
  * the GUI kept a single file per method, so the second species replaced the first;
  * the conservancy compared every species against the same protein set.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from poa.core.parsers import bepipred, netctl, others, papimed
from poa.gui import backend
from poa.services import conservancy_client as cc


# --------------------------------------------------------------------------- header convention
def test_bp3_header_is_protein_specie_id(tmp_path):
    f = tmp_path / "bp3.fasta"
    f.write_text(">E_DENV1_ref\nmktAYIakqr\n")
    df = bepipred.bp3_FastaAnalysis(str(f))
    assert df.iloc[0]["Protein"] == "E"
    assert df.iloc[0]["Specie"] == "DENV1"


def test_bp3_rejects_malformed_header(tmp_path):
    f = tmp_path / "bp3.fasta"
    f.write_text(">Sequence\nmktAYIakqr\n")
    with pytest.raises(ValueError, match="Protein_Specie_ID"):
        bepipred.bp3_FastaAnalysis(str(f))


# --------------------------------------------------------------------------- merging helpers
def _b2(seq, pred):
    return json.dumps({"antigens": {"Sequence": {"AA": list(seq), "PRED": pred}}}).encode("utf-8")


def test_merge_bepipred2_keeps_every_species(tmp_path):
    ctx = backend.WorkContext(tmp_path)
    pred = [0.1, 0.9, 0.9, 0.9, 0.1, 0.1]
    a = backend.import_bepipred2(ctx, "bepipred_denv1.json", _b2("MKTAYI", pred), "DENV1", "E")
    b = backend.import_bepipred2(ctx, "bepipred_chikv.json", _b2("MQRSTA", pred), "CHIKV", "E1")

    out = backend.merge_bepipred2([a.b2_json, b.b2_json], str(tmp_path / "merged.json"))
    antigens = json.loads(Path(out).read_text())["antigens"]
    assert set(antigens) == {"E_DENV1_ref", "E1_CHIKV_ref"}


def test_merge_fastas_dedupes_headers(tmp_path):
    f1 = tmp_path / "a.fasta"
    f1.write_text(">E_DENV1_ref\nMKTA\n")
    f2 = tmp_path / "b.fasta"
    f2.write_text(">E1_CHIKV_ref\nMQRS\n>E_DENV1_ref\nMKTA\n")

    out = backend.merge_fastas([str(f1), str(f2)], str(tmp_path / "all.fasta"))
    headers = [line for line in Path(out).read_text().splitlines() if line.startswith(">")]
    assert headers == [">E_DENV1_ref", ">E1_CHIKV_ref"]


def test_merge_text_concatenates_and_pads_newline(tmp_path):
    f1 = tmp_path / "a.txt"
    f1.write_text(">E_DENV1_ref")          # no trailing newline
    f2 = tmp_path / "b.txt"
    f2.write_text(">E1_CHIKV_ref\n")

    out = backend.merge_text([str(f1), str(f2)], str(tmp_path / "m.txt"))
    assert Path(out).read_text().splitlines() == [">E_DENV1_ref", ">E1_CHIKV_ref"]


def test_consolidate_sources_merges_lists_and_passes_dirs(tmp_path):
    f1 = tmp_path / "a.fasta"
    f1.write_text(">E_DENV1_ref\nMKTA\n")
    f2 = tmp_path / "b.fasta"
    f2.write_text(">E1_CHIKV_ref\nMQRS\n")

    files = backend.consolidate_sources(
        {"b3": [str(f1), str(f2)], "x": [str(f1)], "m": [str(tmp_path)], "p": []}, tmp_path)

    assert "p" not in files                       # empty lists are dropped
    assert files["m"] == str(tmp_path)            # MHC-II directory passes through untouched
    assert files["x"] == str(f1)                  # a single file is used as-is
    assert Path(files["b3"]).name == "bepipred3_merged.fasta"
    assert Path(files["b3"]).read_text().count(">") == 2


# --------------------------------------------------------------------------- species diagnostics
def test_species_diagnostics_flags_missing_reference(tmp_path):
    ref = tmp_path / "ref.fasta"
    ref.write_text(">E_DENV1_ref\nMKTA\n")
    df = pd.DataFrame({"Specie": ["DENV1", "CHIKV"], "Peptide Sequence": ["MK", "TA"]})

    diag = backend.species_diagnostics(str(ref), df)
    assert diag["reference"] == ["DENV1"]
    assert diag["predicted"] == ["CHIKV", "DENV1"]
    assert diag["missing_from_reference"] == ["CHIKV"]
    assert diag["unpredicted"] == []


# --------------------------------------------------------------------------- per-species conservancy
def test_conservancy_uses_each_species_own_proteins(tmp_path):
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    poa1_dir.mkdir(parents=True)
    (poa1_dir / "DENV1_epitopes.fasta").write_text(">DENV1_E_Bepipred2.0_1_3\nTAY\n")
    (poa1_dir / "CHIKV_epitopes.fasta").write_text(">CHIKV_E1_Bepipred2.0_1_3\nRST\n")

    proteins = [("E_DENV1_ref", "MKTAYIA"), ("E1_CHIKV_ref", "MQRSTAY")]
    out_dir = tmp_path / "csvs"
    written = cc.run_conservancy_for_dir(str(poa1_dir), proteins, 70, str(out_dir))
    assert len(written) == 2

    # each epitope is compared against its own species only -> 1 protein in the denominator
    assert "100.00% (1/1)" in (out_dir / "DENV1_conservancy.csv").read_text()
    assert "100.00% (1/1)" in (out_dir / "CHIKV_conservancy.csv").read_text()


def test_conservancy_explicit_specie_set_wins(tmp_path):
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    poa1_dir.mkdir(parents=True)
    (poa1_dir / "DENV1_epitopes.fasta").write_text(">DENV1_E_Bepipred2.0_1_3\nTAY\n")

    world = tmp_path / "denv1_world.fasta"
    world.write_text(">acc1\nMKTAYIA\n>acc2\nMKTGGIA\n")   # identities 100% and 33%
    out_dir = tmp_path / "csvs"

    cc.run_conservancy_for_dir(str(poa1_dir), [("E_DENV1_ref", "MKTAYIA")], 70, str(out_dir),
                               specie_proteins={"DENV1": str(world)})
    assert "50.00% (1/2)" in (out_dir / "DENV1_conservancy.csv").read_text()


def test_conservancy_falls_back_to_full_set_without_species(tmp_path):
    """The realdata CLI passes a world FASTA with accession headers (no species token)."""
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    poa1_dir.mkdir(parents=True)
    (poa1_dir / "DENV1_epitopes.fasta").write_text(">DENV1_E_Bepipred2.0_1_3\nTAY\n")

    out_dir = tmp_path / "csvs"
    cc.run_conservancy_for_dir(str(poa1_dir), [("acc1", "MKTAYIA"), ("acc2", "MKTGGIA")],
                               70, str(out_dir))
    assert "50.00% (1/2)" in (out_dir / "DENV1_conservancy.csv").read_text()


# --------------------------------------------------------------------------- end to end
def test_two_species_in_one_poa1_run(tmp_path):
    ctx = backend.WorkContext(tmp_path / "work")
    pred = [0.1, 0.2, 0.9, 0.9, 0.9, 0.2, 0.1, 0.1]
    prepared = [
        backend.import_bepipred2(ctx, "bepipred_denv1.json", _b2("MKTAYIAK", pred), "DENV1", "E"),
        backend.import_bepipred2(ctx, "bepipred_chikv.json", _b2("MQRSTAYW", pred), "CHIKV", "E1"),
    ]
    proteins = backend.merge_fastas([p.reference_fasta for p in prepared],
                                    str(ctx.inputs_dir / "proteins_all.fasta"))
    files = backend.consolidate_sources({"b2": [p.b2_json for p in prepared]}, ctx.inputs_dir)

    result = backend.run_poa1(ctx, files, {}, proteins)

    assert sorted(result.predictions["Specie"].unique()) == ["CHIKV", "DENV1"]
    assert sorted(result.predictions["Protein"].unique()) == ["E", "E1"]
    assert (ctx.conservancy_epitopes_dir / "DENV1_epitopes.fasta").exists()
    assert (ctx.conservancy_epitopes_dir / "CHIKV_epitopes.fasta").exists()

    report = Path(result.report_path).read_text(encoding="utf-8")
    assert "#Species present in the analysis:\n2" in report

    written = cc.run_conservancy_for_dir(str(ctx.conservancy_epitopes_dir), proteins, 70,
                                         str(ctx.conservancy_csv_dir))
    assert len(written) == 2


# --------------------------------------------------------------------------- one file per species
def test_netctl_two_species_merge(tmp_path):
    a = tmp_path / "netctl_denv1.html"
    a.write_text("NetCTL 1.2\n1 id E_DENV1 pep AADEFGHIK aff 0.1 0.2 cle 0.3 tap 0.4 t2 COMB 0.9 <-E\n")
    b = tmp_path / "netctl_chikv.html"
    b.write_text("NetCTL 1.2\n1 id E1_CHIKV pep BBCDEFGHI aff 0.1 0.2 cle 0.3 tap 0.4 t2 COMB 0.9 <-E\n")

    merged = backend.merge_text([str(a), str(b)], str(tmp_path / "netctl.html"))
    df = netctl.netctlAntigenEpitopes(merged)
    assert sorted(df["Specie"]) == ["CHIKV", "DENV1"]
    assert sorted(df["Protein"]) == ["E", "E1"]


def test_papimed_two_species_merge(tmp_path):
    a = tmp_path / "pap_denv1.txt"
    a.write_text(">E_DENV1_ref\nn\tStart\tSeq\tEnd\n1\t3\tTAYIA\t7\n")
    b = tmp_path / "pap_chikv.txt"
    b.write_text(">E1_CHIKV_ref\nn\tStart\tSeq\tEnd\n1\t3\tRSTAY\t7\n")

    merged = backend.merge_text([str(a), str(b)], str(tmp_path / "pap.txt"))
    df = papimed.PAPepitopes(merged, 0, 0)
    assert sorted(df["Specie"]) == ["CHIKV", "DENV1"]


def test_others_two_species_merge(tmp_path):
    a = tmp_path / "x_denv1.fasta"
    a.write_text(">E_DENV1_OTHER_ref_3_7\nTAYIA\n")
    b = tmp_path / "x_chikv.fasta"
    b.write_text(">E1_CHIKV_OTHER_ref_3_7\nRSTAY\n")

    merged = backend.merge_text([str(a), str(b)], str(tmp_path / "x.fasta"))
    df = others.fasta_epitopes(merged, 0, 0)
    assert sorted(df["Specie"]) == ["CHIKV", "DENV1"]


def test_two_species_inside_one_bepipred2_json(tmp_path):
    """A single BepiPred run over several organisms: keys like 'denv1'/'denv2' carry no protein
    and no ID, so the Protein_Specie_ID regex drops both species unless each is mapped."""
    ctx = backend.WorkContext(tmp_path / "work")
    pred = [0.1, 0.2, 0.9, 0.9, 0.9, 0.2, 0.1, 0.1]
    raw = json.dumps({"antigens": {
        "denv1": {"AA": list("MKTAYIAK"), "PRED": pred},
        "denv2": {"AA": list("MQRSTAYW"), "PRED": pred},
    }}).encode("utf-8")

    raw_path = backend.save_upload(ctx.inputs_dir, "dv1_dv2.json", raw)
    assert backend.antigen_keys(raw_path) == ["denv1", "denv2"]

    mapping = {"denv1": ("NS1", "DENV1", "ref"), "denv2": ("NS1", "DENV2", "ref")}
    prep = backend.import_bepipred2(ctx, "dv1_dv2.json", raw, specie="DENV1", protein="NS1",
                                    mapping=mapping)
    assert [e.header for e in prep.entries] == ["NS1_DENV1_ref", "NS1_DENV2_ref"]
    assert prep.species == ["DENV1", "DENV2"]

    files = backend.consolidate_sources({"b2": [prep.b2_json]}, ctx.inputs_dir)
    result = backend.run_poa1(ctx, files, {}, prep.reference_fasta)
    assert sorted(result.predictions["Specie"].unique()) == ["DENV1", "DENV2"]
    assert list(result.predictions["Protein"].unique()) == ["NS1"]


def test_prepare_bepipred2_without_mapping_is_unchanged(tmp_path):
    """The single-species path (one JSON per virus, used by the realdata CLI) keeps its naming."""
    ctx = backend.WorkContext(tmp_path / "work")
    prep = backend.import_bepipred2(ctx, "bepipred_denv1.json",
                                    _b2("MKTAYI", [0.1, 0.9, 0.9, 0.9, 0.1, 0.1]), "DENV1", "E")
    assert prep.header == "E_DENV1_ref"
    assert Path(prep.b2_json).name == "DENV1_b2.json"
    assert Path(prep.reference_fasta).name == "DENV1_ref.fasta"
    assert prep.species == ["DENV1"]


def test_nonconforming_antigen_keys_gate(tmp_path):
    """The GUI checks the JSON at upload time instead of accepting it and yielding NaN species."""
    ctx = backend.WorkContext(tmp_path / "work")
    pred = [0.1, 0.9, 0.9, 0.9, 0.1, 0.1]

    bad = backend.save_upload(ctx.inputs_dir, "bad.json", json.dumps({"antigens": {
        "denv1": {"AA": list("MKTAYI"), "PRED": pred},
        "Sequence": {"AA": list("MQRSTA"), "PRED": pred},
    }}).encode("utf-8"))
    assert backend.nonconforming_antigen_keys(bad) == ["denv1", "Sequence"]

    good = backend.save_upload(ctx.inputs_dir, "good.json", json.dumps({"antigens": {
        "NS1_DENV1_ref": {"AA": list("MKTAYI"), "PRED": pred},
        "NS1_DENV2_ref": {"AA": list("MQRSTA"), "PRED": pred},
    }}).encode("utf-8"))
    assert backend.nonconforming_antigen_keys(good) == []

    # a conforming file resolves both species without the adapter — this is the path the GUI
    # takes when it accepts the upload directly
    ref = ctx.inputs_dir / "ref.fasta"
    ref.write_text(">NS1_DENV1_ref\nMKTAYI\n>NS1_DENV2_ref\nMQRSTA\n")
    result = backend.run_poa1(ctx, {"b2": good}, {}, str(ref))
    assert sorted(result.predictions["Specie"].unique()) == ["DENV1", "DENV2"]


# --------------------------------------------------------------------------- header inspection
def test_parse_fasta_headers_shows_how_fields_are_read(tmp_path):
    f = tmp_path / "ref.fasta"
    f.write_text(">denv1_ns1\nMKTA\n>Sequence\nMQRS\n")
    rows = backend.parse_fasta_headers(str(f))
    # 'denv1_ns1' is Specie_Protein, so it is read the wrong way round
    assert rows[0]["protein"] == "DENV1" and rows[0]["specie"] == "NS1" and rows[0]["ok"]
    # no separator at all -> no species can be extracted
    assert rows[1]["specie"] == "" and rows[1]["ok"] is False


def test_normalize_fasta_headers_swaps_and_adds_id(tmp_path):
    f = tmp_path / "ref.fasta"
    f.write_text(">denv1_ns1\nMKTA\n>denv2_ns1\nMQRS\n")
    out = backend.normalize_fasta_headers([str(f)], str(tmp_path / "norm.fasta"), swap=True)
    assert [line for line in Path(out).read_text().splitlines() if line.startswith(">")] == [
        ">NS1_DENV1_ref", ">NS1_DENV2_ref"]
    assert backend.species_in_fasta(out) == ["DENV1", "DENV2"]
    # sequences are untouched
    assert "MKTA" in Path(out).read_text() and "MQRS" in Path(out).read_text()


def test_bepipred3_two_species_merge(tmp_path):
    a = tmp_path / "bp3_denv1.fasta"
    a.write_text(">E_DENV1_ref\nmktAYIakq\n")
    b = tmp_path / "bp3_chikv.fasta"
    b.write_text(">E1_CHIKV_ref\nmqrSTAyww\n")

    merged = backend.consolidate_sources({"b3": [str(a), str(b)]}, tmp_path)["b3"]
    df = bepipred.bp3_FastaAnalysis(merged)
    assert sorted(df["Specie"]) == ["CHIKV", "DENV1"]
    assert sorted(df["Protein"]) == ["E", "E1"]
