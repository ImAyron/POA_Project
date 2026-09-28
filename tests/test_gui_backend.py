"""Tests for the GUI backend (UI-agnostic orchestration)."""
import json
from pathlib import Path

import pandas as pd

from poa.gui import backend, viz


def test_default_results_root_is_inside_the_project(monkeypatch, tmp_path):
    monkeypatch.delenv("POA_RESULTS_DIR", raising=False)
    root = backend.default_results_root()
    assert root == backend.REPO_ROOT / "results"
    assert (backend.REPO_ROOT / "poa" / "gui" / "backend.py").exists()  # REPO_ROOT is the repo

    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path / "elsewhere"))
    assert backend.default_results_root() == tmp_path / "elsewhere"


def test_work_context_creates_the_whole_tree(tmp_path):
    ctx = backend.WorkContext(tmp_path / "results")
    for d in (ctx.inputs_dir, ctx.mhcii_dir, ctx.poa1_out_dir,
              ctx.conservancy_csv_dir, ctx.poa2_out_dir, ctx.cache_dir):
        assert d.is_dir()
        assert str(d).startswith(str(ctx.root))


def test_resumable_run_finds_previous_outputs(tmp_path):
    ctx = backend.WorkContext(tmp_path / "results")
    assert backend.resumable_run(ctx) == {"csvs": [], "references": []}

    (ctx.conservancy_csv_dir / "DENV1_conservancy.csv").write_text("Epitope #\n")
    (ctx.conservancy_csv_dir / "notes.txt").write_text("ignored")
    (ctx.inputs_dir / "proteins_all.fasta").write_text(">NS1_DENV1_ref\nMKTA\n")

    found = backend.resumable_run(ctx)
    assert [p.name for p in found["csvs"]] == ["DENV1_conservancy.csv"]
    assert [p.name for p in found["references"]] == ["proteins_all.fasta"]


def test_poa1_result_rebuilt_from_epitope_fastas(tmp_path):
    ctx = backend.WorkContext(tmp_path / "results")
    assert backend.poa1_result_from_disk(ctx) is None      # nothing on disk yet

    ctx.conservancy_epitopes_dir.mkdir(parents=True, exist_ok=True)
    (ctx.conservancy_epitopes_dir / "DENV1_epitopes.fasta").write_text(
        ">DENV1_NS1_Bepipred2.0_3_5\nTAY\n>DENV1_NS1_Bepipred2.0_9_11\nQRQ\n")
    (ctx.conservancy_epitopes_dir / "DENV2_epitopes.fasta").write_text(
        ">DENV2_NS1_Bepipred2.0_3_5\nRST\n")

    result = backend.poa1_result_from_disk(ctx)
    df = result.predictions
    assert len(df) == 3
    assert sorted(df["Specie"].unique()) == ["DENV1", "DENV2"]
    assert list(df["Protein"].unique()) == ["NS1"]
    row = df.iloc[0]
    assert row["Method"] == "Bepipred2.0"
    assert (row["Initial Position"], row["Final Position"]) == ("3", "5")
    assert row["Peptide Sequence"] == "TAY"
    # the rebuilt frame feeds the 2D/3D views directly
    assert viz.epitope_ranges(df, protein="NS1", specie="DENV1") == [(3, 5), (9, 11)]


def test_build_poa1_args_defaults(tmp_path):
    args = backend.build_poa1_args(
        files={"b3": "bp3.fasta"},
        params={"bmin": 5, "mhla": "DP", "e": "y"},
        proteins_path="prot.fasta",
        out_dir=str(tmp_path),
    )
    assert args.b3 == "bp3.fasta"
    assert args.b2 == ""            # unset methods default to ''
    assert args.bmin == 5
    assert args.mhla == "DP"
    assert args.e == "y"
    assert args.f == "prot.fasta"
    assert args.d == str(tmp_path)
    # all attributes the pipeline reads must exist
    for attr in ["b2", "b3", "bmin", "bmax", "p", "pmin", "pmax", "n", "m",
                 "mhla", "mic", "x", "xmin", "xmax", "d", "f", "e"]:
        assert hasattr(args, attr)


def test_build_poa2_args_objective_mapping(tmp_path):
    conserved = backend.build_poa2_args({"objective": "conserved", "t": 70, "rf": 1},
                                        str(tmp_path), "prot.fasta", str(tmp_path))
    assert conserved.g == "True" and conserved.l is None
    assert conserved.t == 70 and conserved.rf == 1

    unique = backend.build_poa2_args({"objective": "unique", "t": 80, "rf": ""},
                                     str(tmp_path), "prot.fasta", str(tmp_path))
    assert unique.g is None and unique.l == "True"
    assert unique.rf is None  # empty -> None


def test_build_poa2_args_threshold_switches(tmp_path):
    default = backend.build_poa2_args({"objective": "conserved", "t": 70},
                                      str(tmp_path), "prot.fasta", str(tmp_path))
    assert default.idf is False and default.strict is False   # historical behaviour by default

    opted_in = backend.build_poa2_args({"objective": "conserved", "t": 70, "idf": True, "strict": True},
                                       str(tmp_path), "prot.fasta", str(tmp_path))
    assert opted_in.idf is True and opted_in.strict is True


def test_count_proteins_by_specie(tmp_path):
    """species_in_fasta de-duplicates; the size of a comparison set is what the threshold needs."""
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(">NS1_DENV1_a\nAAA\n>E_DENV1_b\nCCC\n>E_DENV2_c\nGGG\n")

    assert backend.species_in_fasta(str(fasta)) == ["DENV1", "DENV2"]
    assert backend.count_proteins_by_specie(str(fasta)) == {"DENV1": 2, "DENV2": 1}


def test_import_bepipred2_adapts_upload(tmp_path):
    ctx = backend.WorkContext(tmp_path)
    raw = json.dumps({
        "info": {"failedjobs": 0, "size": 1},
        "antigens": {"Sequence": {"AA": list("MRCVGIGN"), "PRED": [0.1, 0.9, 0.9, 0.9, 0.1, 0.1, 0.9, 0.9]}},
    }).encode("utf-8")

    prep = backend.import_bepipred2(ctx, "bepipred_denv1.json", raw, "DENV1", protein="E")

    assert prep.header == "E_DENV1_ref"
    assert prep.sequence == "MRCVGIGN"
    assert prep.has_x is False
    # both the rewritten JSON and the rebuilt reference land under the session inputs dir
    assert prep.b2_json.startswith(str(ctx.inputs_dir))
    assert prep.reference_fasta.startswith(str(ctx.inputs_dir))
    ref = Path(prep.reference_fasta).read_text().splitlines()
    assert ref[0] == ">E_DENV1_ref" and ref[1] == "MRCVGIGN"


def test_predictions_summary():
    df = pd.DataFrame({
        "Method": ["Bepipred3.0", "Bepipred3.0", "NetCTL"],
        "Specie": ["SARS", "SARS", "MERS"],
        "Peptide Sequence": ["AYI", "LMN", "AADEFGHIK"],
    })
    summ = backend.predictions_summary(df)
    by_method = dict(zip(summ["by_method"]["Method"], summ["by_method"]["count"]))
    assert by_method == {"Bepipred3.0": 2, "NetCTL": 1}
    by_species = dict(zip(summ["by_species"]["Specie"], summ["by_species"]["count"]))
    assert by_species == {"SARS": 2, "MERS": 1}


def test_predictions_summary_empty():
    summ = backend.predictions_summary(pd.DataFrame())
    assert summ["by_method"].empty


def test_run_poa1_via_context(tmp_path):
    ctx = backend.WorkContext(tmp_path / "work")
    proteins = ctx.inputs_dir / "proteins.fasta"
    proteins.write_text(">SPIKE_SARS_NP1\nMKTAYIAMKGVLMNKQRST\n")
    b3 = ctx.inputs_dir / "bp3.fasta"
    b3.write_text(">SPIKE_SARS_NP1\nmktAYIamkgvLMNkqrst\n")

    result = backend.run_poa1(ctx, {"b3": str(b3)}, {}, str(proteins))
    assert len(result.predictions) == 2
    assert (ctx.conservancy_epitopes_dir / "SARS_epitopes.fasta").exists()
