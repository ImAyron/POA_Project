"""Tests for the per-step output previews added to the GUI.

Each POA step now shows what it produced instead of only confirming that it ran. Two of those
previews are not merely cosmetic and are pinned here:

* the step-5 funnel reapplies POA2's own conservancy filters, so its last stage must agree with
  what ``EpitConservAnalysis`` actually selects — otherwise the GUI would explain a selection
  that never happened;
* the step-2 preview deliberately does *not* run the parsers, so it must report what the file
  carries (records, species) and flag headers the pipeline cannot read.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from poa.core.parsers import conservancy as cons
from poa.gui import backend


# --------------------------------------------------------------------------- fixtures
def _csv(directory, specie, rows, threshold=70, operator=">="):
    """Write one IEDB-shaped conservancy CSV; ``rows`` is (name, seq, percent, min, max)."""
    percent_col = cons.percent_column(operator, threshold)
    frame = pd.DataFrame([{
        "Epitope #": i + 1,
        "Epitope name": name,
        "Epitope sequence": seq,
        "Epitope length": len(seq),
        percent_col: percent,
        "Minimum identity": minimum,
        "Maximum identity": maximum,
        "View details": "-",
    } for i, (name, seq, percent, minimum, maximum) in enumerate(rows)])
    path = directory / f"{specie}_conservancy.csv"
    frame.to_csv(path, index=False)
    return path


@pytest.fixture
def csv_dir(tmp_path):
    """Two species whose epitopes land on different sides of the usual filters."""
    directory = tmp_path / "conservancy_csv"
    directory.mkdir()
    _csv(directory, "DENV1", [
        ("DENV1_E_Bepipred2.0_6_93", "INWKGRELK", "100.00% (1/1)", "80.68%", "80.68%"),
        ("DENV1_E_Bepipred2.0_109_120", "MEHKYSWKSWGK", "100.00% (1/1)", "83.33%", "83.33%"),
        ("DENV1_E_Bepipred2.0_124_136", "IGADIQNTT", "0.00% (0/1)", "53.85%", "53.85%"),
    ])
    _csv(directory, "DENV2", [
        ("DENV2_E_Bepipred2.0_140_156", "TAECPNTNR", "100.00% (1/1)", "76.47%", "76.47%"),
        ("DENV2_E_Bepipred2.0_6_136", "VSWKNKELK", "0.00% (0/1)", "67.94%", "67.94%"),
    ])
    return directory


# --------------------------------------------------------------------------- step 5 funnel
@pytest.mark.parametrize("params", [
    dict(objective="conserved", t=70, imin=0, imax=100, m=0, idf=False),
    dict(objective="conserved", t=70, imin=60, imax=100, m=60, idf=False),
    dict(objective="conserved", t=70, imin=70, imax=90, m=100, idf=False),
    dict(objective="conserved", t=70, imin=0, imax=100, m=0, idf=True),
    dict(objective="unique", t=70, imin=0, imax=100, m=0, idf=True),
    dict(objective="conserved", t=70, imin=95, imax=100, m=0, idf=False),  # selects nothing
])
def test_funnel_last_stage_equals_what_poa2_selects(csv_dir, params):
    """The funnel mirrors POA2's filter chain; a divergence between the two fails here."""
    funnel = backend.conservancy_funnel(csv_dir, params)
    symbol = ">=" if params["objective"] == "conserved" else "<"
    selected = cons.EpitConservAnalysis(
        params["t"], symbol, params["m"], params["imax"], params["imin"], str(csv_dir),
        identity_filter=params["idf"])

    assert int(funnel["Epítopos"].iloc[-1]) == len(selected)


def test_funnel_reports_every_filter_in_order(csv_dir):
    funnel = backend.conservancy_funnel(
        csv_dir, dict(objective="conserved", t=70, imin=60, imax=100, m=60, idf=False))

    assert list(funnel.columns) == ["Etapa", "Epítopos"]
    assert int(funnel["Epítopos"].iloc[0]) == 5          # every epitope in both CSVs
    assert funnel["Epítopos"].is_monotonic_decreasing    # a filter can only remove
    # imin=60 drops the 53.85% epitope; m=60 then drops the one 0.00% match left among the four
    assert list(funnel["Epítopos"]) == [5, 4, 4, 3]


def test_funnel_names_the_threshold_it_applied(csv_dir):
    funnel = backend.conservancy_funnel(
        csv_dir, dict(objective="unique", t=70, imin=10, imax=90, m=50, idf=True))
    etapas = list(funnel["Etapa"])
    assert "Identidade mínima ≥ 10%" in etapas
    assert "Identidade máxima ≤ 90%" in etapas
    assert "Proteínas com match ≥ 50%" in etapas
    assert "Identidade do epítopo < 70%" in etapas


def test_funnel_is_empty_without_csvs(tmp_path):
    assert backend.conservancy_funnel(tmp_path, dict(t=70)).empty


def test_previews_skip_csvs_whose_header_cannot_be_read(tmp_path):
    """
    A legacy CSV with no 'at identity ...' header must not break steps 4 and 5.

    POA2 warns about such files and still runs; a preview that raised would turn a warning into
    a dead step, which is how this was first caught.
    """
    directory = tmp_path / "conservancy_csv"
    directory.mkdir()
    pd.DataFrame([{
        "Epitope #": 1, "Epitope name": "DENV1_E_1_3", "Epitope sequence": "AYI",
        "Epitope length": 3, "Percent of protein sequence matches": "100.00% (1/1)",
        "Minimum identity": "80.00%", "Maximum identity": "80.00%", "View details": "-",
    }]).to_csv(directory / "DENV1_conservancy.csv", index=False)

    assert cons.stage_counts(str(directory)) == []
    assert backend.conservancy_funnel(directory, dict(t=70)).empty
    assert backend.conservancy_overview(directory).empty


def test_stage_counts_does_not_raise_on_threshold_mismatch(csv_dir):
    """It describes the data; refusing a mismatch is EpitConservAnalysis's job, not the preview's."""
    stages = cons.stage_counts(str(csv_dir), ID_threshold=95, symbol="<")
    assert stages[0]["remaining"] == 5


# --------------------------------------------------------------------------- step 4 overview
def test_conservancy_overview_counts_epitopes_and_matches(csv_dir):
    overview = backend.conservancy_overview(csv_dir).set_index("Espécie")

    assert overview.loc["DENV1", "Epítopos"] == 3
    assert overview.loc["DENV1", "Com match"] == 2       # the 0.00% row has none
    assert overview.loc["DENV2", "Epítopos"] == 2
    assert overview.loc["DENV2", "Com match"] == 1
    assert overview.loc["DENV1", "Identidade mín. (%)"] == 53.85
    assert overview.loc["DENV1", "Identidade máx. (%)"] == 83.33


def test_conservancy_overview_is_empty_without_csvs(tmp_path):
    assert backend.conservancy_overview(tmp_path).empty


# --------------------------------------------------------------------------- step 2 previews
def test_fasta_preview_reports_records_species_and_proteins(tmp_path):
    path = tmp_path / "bepipred3.fasta"
    path.write_text(">NS1_DENV1_001\nMKTAYI\n>E_DENV2_002\nGGAYIGG\n", encoding="utf-8")

    prev = backend.summarize_source(str(path), "b3")

    assert prev.label == "BepiPred-3.0"
    assert prev.kind == "fasta"
    assert prev.n_records == 2
    assert prev.species == ["DENV1", "DENV2"]
    assert prev.proteins == ["E", "NS1"]
    assert prev.rows[0] == {"Cabeçalho": "NS1_DENV1_001", "Proteína": "NS1",
                            "Espécie": "DENV1", "Resíduos": 6}
    assert prev.notes == []


def test_fasta_preview_flags_headers_with_no_species(tmp_path):
    """The failure the preview exists to catch: POA1 would group these without a species."""
    path = tmp_path / "others.fasta"
    path.write_text(">sequence1\nMKTAYI\n", encoding="utf-8")

    prev = backend.summarize_source(str(path), "x")

    assert prev.species == []
    assert any("sem `_`" in note for note in prev.notes)


def test_json_preview_flags_nonconforming_antigen_keys(tmp_path):
    path = tmp_path / "bp2.json"
    path.write_text(json.dumps({"antigens": {"E_DENV1_ref": {}, "denv2": {}}}), encoding="utf-8")

    prev = backend.summarize_source(str(path), "b2")

    assert prev.kind == "json"
    assert prev.n_records == 2
    assert any("fora da convenção" in note for note in prev.notes)
    conventions = {row["Antígeno"]: row["Convenção"] for row in prev.rows}
    assert conventions == {"E_DENV1_ref": "✅", "denv2": "❌"}


def test_text_preview_shows_the_first_lines(tmp_path):
    path = tmp_path / "pap.txt"
    path.write_text("\n".join(f"line {i}" for i in range(40)), encoding="utf-8")

    prev = backend.summarize_source(str(path), "p")

    assert prev.n_records == 40
    assert prev.head.splitlines() == [f"line {i}" for i in range(15)]


def test_directory_preview_lists_the_mhcii_files(tmp_path):
    directory = tmp_path / "mhcii"
    directory.mkdir()
    (directory / "E_DENV1.html").write_text("<html></html>", encoding="utf-8")
    (directory / "notes.txt").write_text("ignored", encoding="utf-8")

    prev = backend.summarize_source(str(directory), "m")

    assert prev.kind == "dir"
    assert prev.n_records == 1
    assert prev.rows[0]["Arquivo"] == "E_DENV1.html"


def test_missing_file_is_reported_not_raised(tmp_path):
    prev = backend.summarize_source(str(tmp_path / "gone.fasta"), "b3")
    assert prev.notes == ["Arquivo não encontrado."]


def test_summarize_sources_covers_every_collected_file(tmp_path):
    a = tmp_path / "b3.fasta"
    a.write_text(">NS1_DENV1_1\nMKT\n", encoding="utf-8")
    b = tmp_path / "x.fasta"
    b.write_text(">E_DENV2_1\nGGA\n", encoding="utf-8")

    previews = backend.summarize_sources({"b3": [str(a)], "x": [str(b)]})

    assert [p.label for p in previews] == ["BepiPred-3.0", "Outros preditores"]


# --------------------------------------------------------------------------- step 1 residues
def test_step5_shows_the_funnel_even_without_pytmhmm(tmp_path, monkeypatch, csv_dir):
    """
    The pyTMHMM gate blocks the run, not the step.

    pyTMHMM has no Windows wheel and this project's documented route is to tune on Windows and
    execute under WSL. The funnel is pure CSV arithmetic, so it has to render on the side where
    POA2 cannot run — otherwise the parameters can only be judged after switching environments.
    """
    apptest = pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")
    from poa.services import tmhmm_client

    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    monkeypatch.setattr(tmhmm_client, "is_available", lambda: False)
    assert csv_dir.parent == tmp_path   # the fixture already wrote them where WorkContext looks

    app = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")
    at = apptest.AppTest.from_file(app, default_timeout=90)
    at.session_state["run_path"] = str(tmp_path)
    at.session_state["conservancy_ready"] = True
    at.session_state["proteins_path"] = "proteins.fasta"
    at.run()
    at.sidebar.radio[0].set_value("5 · POA2 (conservação + topologia)").run()

    assert not at.exception
    assert any("Efeito dos filtros" in str(s.value) for s in at.subheader)
    assert any("pyTMHMM" in str(e.value) for e in at.error)
    assert next(b for b in at.button if "POA2" in b.label).disabled


def test_fasta_headers_report_length_and_ambiguous_residues(tmp_path):
    path = tmp_path / "proteins.fasta"
    path.write_text(">NS1_DENV1_ref\nMKTAYIX\n>E_DENV2_ref\nGGAYIGG\n", encoding="utf-8")

    rows = {r["header"]: r for r in backend.parse_fasta_headers(str(path))}

    assert rows["NS1_DENV1_ref"]["length"] == 7
    assert rows["NS1_DENV1_ref"]["ambiguous"] == "X"   # POA1 refuses epitopes carrying this
    assert rows["E_DENV2_ref"]["ambiguous"] == ""
