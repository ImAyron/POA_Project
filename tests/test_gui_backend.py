"""Tests for the GUI backend (UI-agnostic orchestration)."""
import json
from pathlib import Path

import pandas as pd

from poa.gui import backend


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
    b3.write_text(">SARS_SPIKE_NP1\nmktAYIamkgvLMNkqrst\n")

    result = backend.run_poa1(ctx, {"b3": str(b3)}, {}, str(proteins))
    assert len(result.predictions) == 2
    assert (ctx.conservancy_epitopes_dir / "SARS_epitopes.fasta").exists()
