"""Regression cases found in the repository review; no live predictors required."""
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from poa.core import dbsearch, pipeline, topology
from poa.core.parsers import bepipred, conservancy, others
from poa.gui import backend
from poa.services import database_client as db


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


def test_other_predictor_header_without_ncbi_id_is_accepted(tmp_path):
    """
    The NCBI ID is optional, so five fields is a valid header.

    The README documents it as ``<ID_NCBI(se_houver)>`` — "if there is one" — but the parser
    required six fields and raised on every epitope from a predictor that reports no accession,
    rejecting a whole file the convention allows.
    """
    path = tmp_path / "epitopes.fasta"
    path.write_text(">E_DENV1_OTHER_2_7\nKTAYIA\n", encoding="utf-8")

    result = others.fasta_epitopes(path, 0, 0)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["Protein"] == "E"
    assert row["Specie"] == "DENV1"
    assert row["Method"] == "OTHER"
    assert row["ID_Sequence"] == ""          # absent, not borrowed from the position fields
    assert row["Initial Position"] == "2"
    assert row["Final Position"] == "7"
    assert row["Peptide Sequence"] == "KTAYIA"


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


def test_bepipred2_refuses_an_antigen_key_it_cannot_resolve(tmp_path):
    """
    Regression: a key missing a field resolved to NaN in Protein/Specie/ID_Sequence, and because
    NaN never equals itself the per-species filter matched no row — POA1 reported success with the
    epitope in its table while writing a 0-byte 'nan_epitopes.fasta' for the conservancy stage.
    bp3_FastaAnalysis already refused the same malformed header; the CLI had no guard.
    """
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps({"antigens": {
        "SPIKE_SARS": {"AA": list("MKT"), "PRED": [0.9, 0.9, 0]},
    }}), encoding="utf-8")

    with pytest.raises(ValueError, match="SPIKE_SARS"):
        bepipred.bp2_AntigenEpitopes(bepipred.bp2_JsonAnalysis(path))


def test_topology_disqualifies_a_short_record_instead_of_ending_the_run(tmp_path, monkeypatch):
    """
    Regression: the range check sat inside the candidate loop and raised, so a second record of
    the same protein/species carrying the peptide — a shorter isolate — aborted the whole POA2
    even when the correct longer record had already matched, in either record order.
    """
    path = tmp_path / "proteins.fasta"
    # the epitope declares 4-6, which only the 9-aa record can host
    path.write_text(">E_DENV1_long\nMKTAYIQQQ\n>E_DENV1_short\nAYIQQ\n", encoding="utf-8")
    monkeypatch.setattr(topology, "pyTMHMMpredict", lambda _: (
        ["E_DENV1_LONG", "E_DENV1_SHORT"], ["ooommmooo", "ooooo"]))

    frame = pd.DataFrame({"Epitope name": ["DENV1_E_OTHER_4_6"], "Epitope sequence": ["AYI"]})
    result = topology.tmhmmAnalysis(SimpleNamespace(f=str(path)), frame)

    # the long record is the one that can host the coordinates, so its 'mmm' annotation is used
    assert result.iloc[0]["Portion_TM"] == 1.0
    assert result.iloc[0]["Portion_Outside"] == 0.0


def test_topology_disqualifies_a_short_record_whichever_order_it_appears_in(tmp_path, monkeypatch):
    """The short record first used to abort just the same: the loop had no break before the raise."""
    path = tmp_path / "proteins.fasta"
    path.write_text(">E_DENV1_short\nAYIQQ\n>E_DENV1_long\nMKTAYIQQQ\n", encoding="utf-8")
    monkeypatch.setattr(topology, "pyTMHMMpredict", lambda _: (
        ["E_DENV1_SHORT", "E_DENV1_LONG"], ["ooooo", "ooommmooo"]))

    frame = pd.DataFrame({"Epitope name": ["DENV1_E_OTHER_4_6"], "Epitope sequence": ["AYI"]})
    result = topology.tmhmmAnalysis(SimpleNamespace(f=str(path)), frame)

    assert result.iloc[0]["Portion_TM"] == 1.0


def test_topology_reports_coordinates_that_fit_no_record_without_aborting(tmp_path, monkeypatch):
    """
    Coordinates outside every record carrying the peptide are still a data error worth seeing,
    but they belong to one row: the epitope keeps the '-' placeholder like any other unresolved
    one, and the warning names the positions rather than ending the run.
    """
    path = tmp_path / "proteins.fasta"
    path.write_text(">E_DENV1_short\nAYIQQ\n", encoding="utf-8")
    monkeypatch.setattr(topology, "pyTMHMMpredict", lambda _: (["E_DENV1_SHORT"], ["ooooo"]))

    frame = pd.DataFrame({"Epitope name": ["DENV1_E_OTHER_40_42"], "Epitope sequence": ["AYI"]})
    with pytest.warns(UserWarning, match="fall outside every record"):
        result = topology.tmhmmAnalysis(SimpleNamespace(f=str(path)), frame)

    assert len(result) == 1
    assert result.iloc[0]["Portion_TM"] == "-"
    assert result.iloc[0]["Portion_Outside"] == "-"


def _poa0_args(tmp_path, accession, *, nome, sigla="ZIKA"):
    """The namespace ``run_poa0`` expects, for one UniProt entry."""
    return SimpleNamespace(s="", a=accession, d=str(tmp_path), r="", nome=nome, sigla=sigla,
                           estrutura="", descendentes=False, limite=25, l="")


def test_poa0_keeps_two_isolates_of_one_protein_apart(tmp_path, monkeypatch):
    """
    Two entries of the same protein and species must not land on the same file.

    The header told them apart by accession but the file name did not, so the second fetch
    reopened the first one's path in mode 'w' and replaced its sequence — and because the panel
    deduplicates by path, the screen still showed a single file. Building a diversity set means
    fetching exactly such entries, so the collision struck the one workflow that needs them.
    """
    sequences = {"Q32ZE1": "ACDEFGHIKL" * 10, "A0A024B7W1": "MNPQRSTVWY" * 10}
    monkeypatch.setattr(db, "get_entry", lambda acc: {
        "sequence": {"length": len(sequences[acc]), "value": sequences[acc]},
        "organism": {"scientificName": "Zika virus"},
        "proteinDescription": {"submissionNames": [{"fullName": {"value": "Genome polyprotein"}}]},
        "features": []})

    first = dbsearch.run_poa0(_poa0_args(tmp_path, "Q32ZE1", nome="Genomepolyprotein"))
    second = dbsearch.run_poa0(_poa0_args(tmp_path, "A0A024B7W1", nome="Genomepolyprotein"))

    assert first["arquivos"]["fasta"] != second["arquivos"]["fasta"]
    # Each path still holds the sequence its own fetch produced.
    for report, accession in ((first, "Q32ZE1"), (second, "A0A024B7W1")):
        lines = Path(report["arquivos"]["fasta"]).read_text(encoding="utf-8").splitlines()
        assert lines[0] == f">Genomepolyprotein_ZIKA_{accession}"
        assert "".join(lines[1:]) == sequences[accession]


def test_comparison_set_accumulates_without_counting_a_refetch_twice(tmp_path):
    """
    The POA0 panel sums each fetch into the species' comparison FASTA, in place.

    Reading and rewriting the same path is what a second fetch does, and fetching an accession
    that is already in the set must leave the denominator alone — a doubled protein would inflate
    the conservancy denominator without adding any diversity to it.
    """
    target = tmp_path / "comparison_ZIKA.fasta"
    one = tmp_path / "one.fasta"
    one.write_text(">Genomepolyprotein_ZIKA_Q32ZE1\nMKNPKKKS\n", encoding="utf-8")
    two = tmp_path / "two.fasta"
    two.write_text(">Genomepolyprotein_ZIKA_A0A024B7W1\nMKNPKKKT\n", encoding="utf-8")

    backend.merge_fastas([str(one)], str(target))
    backend.merge_fastas([str(target), str(two)], str(target))
    assert [r["header"] for r in backend.parse_fasta_headers(str(target))] == [
        "Genomepolyprotein_ZIKA_Q32ZE1", "Genomepolyprotein_ZIKA_A0A024B7W1"]

    backend.merge_fastas([str(target), str(two)], str(target))
    assert len(backend.parse_fasta_headers(str(target))) == 2


def test_comparison_set_is_restored_onto_step_4(tmp_path):
    """
    A set the POA0 panel built comes back on the step that consumes it, not as a reference.

    ``comparison_<SPECIE>.fasta`` is the name step 4 already writes for an uploaded set, so the
    provenance recorded here has to classify it the same way — a set returned as a '-f' reference
    would make POA1 predict on every isolate in the diversity set.
    """
    ctx = backend.WorkContext(str(tmp_path))
    ctx.inputs_dir.mkdir(parents=True, exist_ok=True)
    target = ctx.inputs_dir / "comparison_ZIKA.fasta"
    target.write_text(">Genomepolyprotein_ZIKA_Q32ZE1\nMKNPKKKS\n", encoding="utf-8")
    backend.record_input(ctx, target, role="comparison_set", specie="ZIKA",
                         original_name=target.name)

    restored = backend.restore_analysis(ctx)
    assert restored.comparison_sets == {"ZIKA": str(target)}
    assert str(target) not in restored.references
