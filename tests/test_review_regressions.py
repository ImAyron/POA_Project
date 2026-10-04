"""Regression cases found in the repository review; no live predictors required."""
import json
from types import SimpleNamespace

import pandas as pd
import pytest

from poa.core import pipeline, topology
from poa.core.parsers import bepipred, conservancy, others
from poa.gui import backend


def test_bepipred_does_not_join_consecutive_positions_from_different_antigens(tmp_path):
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps({"antigens": {
        "E_DENV1_ref": {"AA": list("MKT"), "PRED": [0, 0.9, 0]},
        "E_DENV2_ref": {"AA": list("GAY"), "PRED": [0, 0, 0.9]},
    }}), encoding="utf-8")
    result = bepipred.bp2_AntigenEpitopes(bepipred.bp2_JsonAnalysis(path))
    assert list(result["Peptide Sequence"]) == ["K", "Y"]
    assert list(result["Specie"]) == ["DENV1", "DENV2"]
    assert list(result["Initial Position"]) == [2, 3]


def test_wrapped_fasta_is_one_epitope_and_keeps_full_id(tmp_path):
    path = tmp_path / "epitopes.fasta"
    path.write_text(">E_DENV1_OTHER_NP_123_2_7 description\nKT\nAYIA\n", encoding="utf-8")
    result = others.fasta_epitopes(path, 6, 6)
    assert len(result) == 1
    assert result.iloc[0]["Peptide Sequence"] == "KTAYIA"
    assert result.iloc[0]["ID_Sequence"] == "NP_123"


def test_invalid_other_predictor_header_has_actionable_error(tmp_path):
    path = tmp_path / "invalid.fasta"
    path.write_text(">sequence\nMKT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Protein_Specie_Method_ID_Init_Final"):
        others.fasta_epitopes(path, 0, 0)


def test_poa1_can_report_no_epitopes_after_filtering(tmp_path):
    path = tmp_path / "proteins.fasta"
    path.write_text(">E_DENV1_ref\nMKTAYIA\n", encoding="utf-8")
    args = backend.build_poa1_args({"b3": str(path)}, {"bmin": 20},
                                   str(path), str(tmp_path / "output"))
    result = pipeline.run_poa1(args)
    assert result.predictions.empty
    report = (tmp_path / "output" / "Analysis_report.txt").read_text(encoding="utf-8")
    assert "Total:\t0" in report
    assert "inf" not in report
    assert "0 aa." in report


def test_topology_uses_named_protein_and_declared_repeat_position(tmp_path, monkeypatch):
    path = tmp_path / "proteins.fasta"
    path.write_text(">NS1_DENV1_ref\nAYIQQQAYI\n>E_DENV1_ref\nAYIQQQAYI\n", encoding="utf-8")
    monkeypatch.setattr(topology, "pyTMHMMpredict", lambda _: (
        ["NS1_DENV1_REF", "E_DENV1_REF"], ["ooooooooo", "oooooommm"]))
    frame = pd.DataFrame({"Epitope name": ["DENV1_E_OTHER_7_9"], "Epitope sequence": ["AYI"]})
    result = topology.tmhmmAnalysis(SimpleNamespace(f=str(path)), frame)
    assert result.iloc[0]["Portion_TM"] == 1.0
    assert result.iloc[0]["Portion_Outside"] == 0.0


def test_empty_topology_does_not_invoke_predictor(tmp_path, monkeypatch):
    path = tmp_path / "proteins.fasta"
    path.write_text(">E_DENV1_ref\nMKT\n", encoding="utf-8")
    def unexpected_predictor(_):
        pytest.fail("No topology prediction is needed for zero epitopes")
    monkeypatch.setattr(topology, "pyTMHMMpredict", unexpected_predictor)
    result = topology.tmhmmAnalysis(SimpleNamespace(f=str(path)),
                                    pd.DataFrame(columns=["Epitope name", "Epitope sequence"]))
    assert result.empty
    assert "Portion_TM" in result.columns


@pytest.mark.parametrize("content", ["", "not,a,conservancy,csv\n1,2,3,4\n",
    "Minimum identity,Percent of protein sequence matches at identity >= 70%\n80%,100%\n"])
def test_bad_csv_does_not_crash_previews_and_reports_filename(tmp_path, content):
    (tmp_path / "bad.csv").write_text(content, encoding="utf-8")
    overview = backend.conservancy_overview(tmp_path)
    assert overview.empty
    assert "bad.csv" in overview.attrs["errors"][0]
    assert backend.conservancy_funnel(tmp_path, {"t": 70}).empty


def test_csv_discovery_is_sorted_and_ignores_backup_files(tmp_path):
    for name in ["b.CSV", "a.csv", "old.csv.bak", "notes.txt"]:
        (tmp_path / name).touch()
    assert conservancy.map_EpConservFiles(tmp_path) == [str(tmp_path / "a.csv"), str(tmp_path / "b.CSV")]


def test_mixed_thresholds_do_not_crash_funnel_but_poa2_rejects_them(tmp_path):
    for threshold in [70, 80]:
        pd.DataFrame([{"Epitope name": "DENV1_E_OTHER_1_3", "Epitope sequence": "MKT",
                       "Epitope length": 3, "Minimum identity": "90%", "Maximum identity": "100%",
                       conservancy.percent_column(">=", threshold): "100% (1/1)"}]).to_csv(
                           tmp_path / f"{threshold}.csv", index=False)
    assert backend.conservancy_funnel(tmp_path, {"t": 70}).empty
    with pytest.raises(ValueError, match="different sequence identity thresholds"):
        conservancy.EpitConservAnalysis(70, ">=", 60, 100, 60, tmp_path)


def test_funnel_defaults_match_execution_defaults(tmp_path):
    pd.DataFrame([{"Epitope name": "DENV1_E_OTHER_1_3", "Epitope sequence": "MKT",
                   "Epitope length": 3, "Minimum identity": "50%", "Maximum identity": "100%",
                   conservancy.percent_column(">=", 70): "100% (1/1)"}]).to_csv(
                       tmp_path / "data.csv", index=False)
    params = {"t": 70}
    args = backend.build_poa2_args(params, str(tmp_path), "unused", "unused")
    selected = conservancy.EpitConservAnalysis(args.t, ">=", args.m, args.imax, args.imin, args.d)
    assert backend.conservancy_funnel(tmp_path, params)["Epítopos"].iloc[-1] == len(selected) == 0
