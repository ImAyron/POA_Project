"""Behavioral regressions for the fifteen findings of the October 15 audit."""
import json
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from poa.cli.poa0 import build_parser
from poa.core import dbsearch, pipeline
from poa.core.artifacts import staged_outputs
from poa.core.models import STANDARD_COLUMNS
from poa.core.parsers import conservancy as cp, mhcii
from poa.gui import backend, viz
from poa.services import conservancy_client as cc, mhcii_client
from tools import validate_conservancy


def put(root, name, text):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def poa1_inputs(root, species="AB"):
    ctx = backend.WorkContext(root)
    refs = put(ctx.inputs_dir, "ref.fasta", ">E_A_ref\nACDE\n>E_B_ref\nACDE\n")
    pred = put(ctx.inputs_dir, "pred.fasta", "".join(f">E_{s}_M_1_4\nACDE\n" for s in species))
    return ctx, refs, pred


def csv_input(root, threshold=70):
    rows = cc.conservancy_rows([("A_E_M_1_4", "ACDE")], [("E_A_ref", "ACDE")], threshold)
    return Path(put(root, "A.csv", cc.rows_to_csv(rows)))


def test_01_rerun_removes_species_and_empty_outputs(tmp_path):
    ctx, refs, pred = poa1_inputs(tmp_path)
    result = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    cc.run_conservancy_for_dir(result.conservancy_dir, refs, 70, str(ctx.conservancy_csv_dir))
    put(ctx.inputs_dir, "pred.fasta", ">E_A_M_1_4\nACDE\n")
    result = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    cc.run_conservancy_for_dir(result.conservancy_dir, refs, 70, str(ctx.conservancy_csv_dir))
    assert not (ctx.conservancy_epitopes_dir / "B_epitopes.fasta").exists()
    assert [Path(p).name for p in cp.map_EpConservFiles(str(ctx.conservancy_csv_dir))] == ["A_conservancy.csv"]
    assert cp.stale_csvs(str(ctx.conservancy_csv_dir)) == []
    assert list(ctx.poa1_out_dir.glob(".poa_history/*/Conservancy Analysis/B_epitopes.fasta"))
    put(ctx.inputs_dir, "pred.fasta", "")
    result = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    assert result.predictions.empty
    assert not list(ctx.conservancy_epitopes_dir.glob("*.fasta"))
    cc.run_conservancy_for_dir(result.conservancy_dir, refs, 70, str(ctx.conservancy_csv_dir))
    assert cp.map_EpConservFiles(str(ctx.conservancy_csv_dir)) == []


def test_failed_poa1_preserves_previous_batch(tmp_path, monkeypatch):
    ctx, refs, pred = poa1_inputs(tmp_path)
    result = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    before = Path(result.report_path).read_bytes()
    monkeypatch.setattr(pipeline.report, "writereport", Mock(side_effect=OSError("disk full")))
    with pytest.raises(OSError, match="disk full"):
        backend.run_poa1(ctx, {"x": pred}, {"e": "y"}, refs)
    assert Path(result.report_path).read_bytes() == before
    assert (ctx.conservancy_epitopes_dir / "B_epitopes.fasta").exists()
    assert not list(ctx.poa1_out_dir.glob("*.xlsx"))


def test_publication_rolls_back_a_failed_replace(tmp_path, monkeypatch):
    from poa.core import artifacts
    put(tmp_path, "a.txt", "old-a")
    put(tmp_path, "b.txt", "old-b")
    put(tmp_path, "user.txt", "keep")
    replace = artifacts.os.replace
    failed = False

    def fail_once(src, dest):
        nonlocal failed
        if not failed and Path(src).parent.name.startswith(".poa_stage_") and Path(src).name == "b.txt":
            failed = True
            raise OSError("publication failed")
        return replace(src, dest)

    monkeypatch.setattr(artifacts.os, "replace", fail_once)
    with pytest.raises(OSError, match="publication failed"):
        with staged_outputs(tmp_path, ["a.txt", "b.txt"]) as stage:
            put(stage, "a.txt", "new-a")
            put(stage, "b.txt", "new-b")
    assert (tmp_path / "a.txt").read_text() == "old-a"
    assert (tmp_path / "b.txt").read_text() == "old-b"
    assert (tmp_path / "user.txt").read_text() == "keep"


def test_02_metadata_cannot_override_csv_header(tmp_path):
    csv_input(tmp_path, 90)
    cp.write_metadata(str(tmp_path), 70, ">=", species={"A": {"csv": "A.csv"}})
    with pytest.raises(ValueError, match="header.*does not match"):
        cp.EpitConservAnalysis(70, ">=", 0, 100, 0, str(tmp_path), strict=True)


def test_metadata_detects_changed_content_with_same_header(tmp_path):
    path = csv_input(tmp_path)
    cp.write_metadata(str(tmp_path), 70, ">=")
    path.write_text(path.read_text().replace("100.00%", "80.00%"))
    with pytest.raises(ValueError, match="changed"):
        cp.inspect_directory(str(tmp_path))


def tsv(sequence):
    return "allele\tseq_num\tstart\tend\tpeptide\tic50\n" + f"HLA-DRB1*01:01\t1\t1\t4\t{sequence}\t10\n"


def test_03_mhcii_references_survive_consolidation(tmp_path):
    ctx, refs, _ = poa1_inputs(tmp_path, "A")
    for ident, sequence in [("ref1", "ACDE"), ("ref2", "AAAA")]:
        mhcii_client.write_result_file(tsv(sequence), "E", "A", str(ctx.mhcii_dir), reference_id=ident)
    result = backend.run_poa1(ctx, {"m": str(ctx.mhcii_dir)}, {}, refs)
    assert set(result.predictions["ID_Sequence"]) == {"ref1", "ref2"}
    assert set(result.predictions["Peptide Sequence"]) == {"ACDE", "AAAA"}


def test_04_conflicting_merge_preserves_previous_destination(tmp_path):
    a = put(tmp_path, "a.fasta", ">E_A_ref\nAAAA\n")
    b = put(tmp_path, "b.fasta", ">E_A_ref\nCCCC\n")
    dest = put(tmp_path, "merged.fasta", "previous output")
    with pytest.raises(ValueError, match="Conflicting FASTA ID"):
        backend.merge_fastas([a, b], dest)
    assert Path(dest).read_text() == "previous output"


def test_05_reference_is_the_poa1_snapshot_not_a_diversity_file(tmp_path):
    ctx, refs, pred = poa1_inputs(tmp_path)
    backend.run_poa1(ctx, {"x": pred}, {}, refs)
    put(ctx.inputs_dir, "comparison_A.fasta", ">E_A_variant\nCCCC\n")
    Path(refs).write_text(">E_A_ref\nAAAA\n")
    assert backend.reference_sequence(ctx, "E", "A", refs)["sequence"] == "ACDE"
    assert backend.reference_sequence(ctx, "E", "A")["path"] == str(ctx.poa1_out_dir / "reference.fasta")


def test_legacy_reference_requires_explicit_path_and_accepts_faa(tmp_path):
    ctx = backend.WorkContext(tmp_path)
    path = put(ctx.inputs_dir, "ref.faa", ">E_A_ref\nACDE\n")
    assert backend.reference_sequence(ctx, "E", "A") is None
    assert backend.reference_sequence(ctx, "E", "A", path)["sequence"] == "ACDE"
    Path(path).write_text(">E_A_ref\nACDE\n>E_A_other\nCCCC\n")
    with pytest.raises(ValueError, match="referências diferentes"):
        backend.reference_sequence(ctx, "E", "A", path)


@pytest.mark.parametrize("region", ["0-3", "4-2", "2-999"])
def test_06_poa0_rejects_out_of_range_regions(tmp_path, monkeypatch, region):
    """
    Uma região fora da sequência para o POA0 antes de ele escrever qualquer coisa.

    O erro é `SystemExit`, e não `ValueError`, porque é erro de entrada do usuário: todo o resto
    de `run_poa0` (acesso irresolvível, busca vazia, sequência ausente) sai por `SystemExit` com
    uma linha em português, e `poa/cli/poa0.py` não tem handler para exceção nenhuma — um
    `ValueError` viraria traceback na cara de quem só errou o `-r`.
    """
    monkeypatch.setattr(dbsearch.db, "get_entry", lambda _: {
        "sequence": {"value": "ACDE", "length": 4}, "organism": {"scientificName": "Species A"}})
    out = tmp_path / "out"
    args = build_parser().parse_args(["-a", "TEST", "-r", region, "--nome", "E", "-d", str(out)])
    with pytest.raises(SystemExit, match="nao cabe em TEST"):
        dbsearch.run_poa0(args)
    assert not out.exists()


def test_07_nested_mhcii_paths_are_readable(tmp_path):
    path = put(tmp_path, "sub/E_A.html", mhcii_client.normalize_tsv(tsv("ACDE")))
    species, proteins, paths = mhcii.files_map_MHCIIBD(str(tmp_path))
    assert paths == [path] and species == ["A"] and proteins == ["E"]
    assert len(mhcii.MHCIIHTMLconverter(paths[0])) == 1


@pytest.mark.parametrize("name", ["../outside.txt", "..\\outside.txt", "/outside.txt", "C:\\outside.txt", "NUL", "file:stream", "file."])
def test_08_upload_paths_are_contained(tmp_path, name):
    with pytest.raises(ValueError, match="filename"):
        backend.save_upload(tmp_path / "uploads", name, b"audit")
    assert not (tmp_path / "outside.txt").exists()


def test_09_create_never_reopens_an_existing_analysis(tmp_path):
    first = backend.create_run(tmp_path, "Teste A")
    put(first.inputs_dir, "keep.txt", "old")
    second = backend.create_run(tmp_path, "Teste-A")
    assert second.root != first.root
    assert (first.inputs_dir / "keep.txt").exists()
    assert not (second.inputs_dir / "keep.txt").exists()
    assert len(backend.list_runs(tmp_path)) == 2


def test_11_viability_uses_matching_species_and_protein(tmp_path):
    ctx, refs, pred = poa1_inputs(tmp_path, "A")
    Path(refs).write_text(">E_A_ref\nACDE\n>E_B_ref\nAC\n>SHORT_A_ref\nA\n")
    assert len(backend.run_poa1(ctx, {"x": pred}, {}, refs).predictions) == 1
    Path(refs).write_text(">E_A_ref\nAC\n>E_B_ref\nACDE\n")
    with pytest.raises(ValueError, match="matching reference"):
        backend.run_poa1(ctx, {"x": pred}, {}, refs)


def test_12_manual_metadata_has_a_valid_inventory(tmp_path):
    csv_input(tmp_path)
    cp.write_metadata(str(tmp_path), 70, ">=", source="manual")
    assert cp.stale_csvs(str(tmp_path)) == []
    assert len(cp.EpitConservAnalysis(70, ">=", 0, 100, 0, str(tmp_path), strict=True)) == 1
    put(tmp_path, "extra.csv", (tmp_path / "A.csv").read_text())
    assert cp.stale_csvs(str(tmp_path)) == ["extra.csv"]


def test_legacy_manual_metadata_without_inventory_is_not_an_empty_inventory(tmp_path):
    csv_input(tmp_path)
    put(tmp_path, cp.METADATA_FILENAME, json.dumps({"threshold": 70, "operator": ">=", "source": "manual"}))
    assert cp.stale_csvs(str(tmp_path)) == []


@pytest.mark.parametrize("number", ["NaN", "inf", "-inf", "", "-1", "101"])
def test_13_validator_refuses_invalid_metrics(tmp_path, number):
    path = csv_input(tmp_path)
    df = pd.read_csv(path, keep_default_na=False)
    df["Minimum identity"] = number
    df.to_csv(path, index=False)
    ep = put(tmp_path, "ep.fasta", ">A_E_M_1_4\nACDE\n")
    ref = put(tmp_path, "ref.fasta", ">E_A_ref\nACDE\n")
    assert validate_conservancy.main(["-i", str(path), "-e", ep, "-p", ref, "-t", "70"]) == 2


def open_app(ctx, refs, pred, monkeypatch, **extra):
    from streamlit.testing.v1 import AppTest
    from poa import logging_conf
    monkeypatch.setattr(logging_conf, "use_log_file", lambda _: None)
    app = Path(__file__).resolve().parents[1] / "poa/gui/app.py"
    at = AppTest.from_file(str(app), default_timeout=90)
    state = {"run_path": str(ctx.root), "proteins_path": refs, "poa1_sources": {"x": [pred]}}
    state.update(extra)
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    assert not at.exception
    return at


def test_14_poa1_rerun_invalidates_poa2_even_on_failure(tmp_path, monkeypatch):
    ctx, refs, pred = poa1_inputs(tmp_path, "A")
    old = pipeline.Poa2Result(pd.DataFrame(), "old.xlsx", None, ">=")
    at = open_app(ctx, refs, pred, monkeypatch, poa2_result=old, conservancy_ready=True)
    at.sidebar.radio[0].set_value("3 · POA1 & Relatório").run()
    next(b for b in at.button if "Executar POA1" in b.label).click().run()
    assert not at.exception
    assert at.session_state["poa2_result"] is None
    assert at.session_state["conservancy_ready"] is False
    at.session_state["poa2_result"] = old
    monkeypatch.setattr(backend, "run_poa1", Mock(side_effect=ValueError("failed")))
    next(b for b in at.button if "Executar POA1" in b.label).click().run()
    assert at.session_state["poa2_result"] is None
    assert at.session_state["poa1_result"] is None


def test_changed_reference_invalidates_results_without_running(tmp_path, monkeypatch):
    ctx, refs, pred = poa1_inputs(tmp_path, "A")
    at = open_app(ctx, refs, pred, monkeypatch)
    at.sidebar.radio[0].set_value("3 · POA1 & Relatório").run()
    next(b for b in at.button if "Executar POA1" in b.label).click().run()
    assert at.session_state["poa1_result"] is not None
    Path(refs).write_text(">E_A_ref\nCCCC\n")
    at.run()
    assert not at.exception
    assert at.session_state["poa1_result"] is None


def pdb(sequence="ACDE", chain="A"):
    names = dict(A="ALA", C="CYS", D="ASP", E="GLU")
    return "\n".join(f"ATOM  {i:>5}  CA  {names[aa]} {chain}{i:>4}    "
                     f"{float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C"
                     for i, aa in enumerate(sequence, 1)) + "\n"


def test_10_viewer_gets_the_validated_chain(tmp_path, monkeypatch):
    ctx, refs, pred = poa1_inputs(tmp_path, "A")
    result = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    path = put(ctx.inputs_dir, "structure.pdb", pdb() + pdb("CCCC", "B"))
    render = Mock(return_value="<div>structure</div>")
    monkeypatch.setattr(viz, "build_3dmol_view_html", render)
    at = open_app(ctx, refs, pred, monkeypatch, poa1_result=result, pdb_path=path)
    at.sidebar.radio[0].set_value("7 · Visualização (2D/3D)").run()
    assert not at.exception
    assert render.call_args.kwargs["chain"] == "A"
    assert render.call_args.args[1] == [(1, 4)]


def test_15_empty_predictions_render_structure_without_highlights(tmp_path, monkeypatch):
    ctx, refs, pred = poa1_inputs(tmp_path, "A")
    result = pipeline.Poa1Result(pd.DataFrame(columns=STANDARD_COLUMNS), "", "")
    path = put(ctx.inputs_dir, "structure.pdb", pdb())
    render = Mock(return_value="<div>structure</div>")
    monkeypatch.setattr(viz, "build_3dmol_view_html", render)
    at = open_app(ctx, refs, pred, monkeypatch, poa1_result=result, pdb_path=path)
    at.sidebar.radio[0].set_value("7 · Visualização (2D/3D)").run()
    assert not at.exception
    assert render.call_args.args[1] == []
