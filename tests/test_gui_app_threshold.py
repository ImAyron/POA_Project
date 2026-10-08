"""The threshold guard as the user meets it, driving the real Streamlit app headlessly.

This logic lives in ``poa/gui/app.py`` and is the user-facing half of the fix: POA2 must take the
threshold from the data rather than from the form, and must not let a run proceed under an
objective the CSVs were not computed with.
"""
import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from pathlib import Path  # noqa: E402

from streamlit.testing.v1 import AppTest  # noqa: E402

from poa.core.parsers import conservancy as cons_parser  # noqa: E402
from poa.services import conservancy_client as cc  # noqa: E402
from poa.services import tmhmm_client  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")
STEP5 = "5 · POA2 (conservação + topologia)"


@pytest.fixture
def results_root(tmp_path, monkeypatch):
    """A results tree holding conservancy CSVs computed at >= 70%, with POA2 reachable."""
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    # POA2 is gated on pyTMHMM, which has no Windows wheel; the guard under test runs before it.
    monkeypatch.setattr(tmhmm_client, "is_available", lambda: True)

    epitopes = tmp_path / "poa1_out" / "Conservancy Analysis"
    epitopes.mkdir(parents=True)   # a pre-run-folder tree: the legacy analysis, used as-is
    (epitopes / "SP_epitopes.fasta").write_text(">SP_E_BEPIPRED_1_5\nAYIAM\n")
    cc.run_conservancy_for_dir(
        str(epitopes),
        [("E_SP_p1", "GGAYIAMGG"), ("E_SP_p2", "TTAYIAMTT"), ("E_SP_p3", "GGAYIALGG")],
        threshold=70, out_dir=str(tmp_path / "conservancy_csv"), operator=">=")
    return tmp_path


def _open_step5(root, **state):
    at = AppTest.from_file(APP, default_timeout=90)
    at.session_state["run_path"] = str(root)   # skip the analysis picker
    at.session_state["conservancy_ready"] = True
    at.session_state["proteins_path"] = "proteins.fasta"
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    at.sidebar.radio[0].set_value(STEP5).run()
    return at


def _run_button(at):
    return next(b for b in at.button if "POA2" in b.label)


def test_threshold_comes_from_the_csvs_not_from_the_form(results_root):
    """The session says 90; the CSVs were built at 70. The data must win."""
    at = _open_step5(results_root, threshold=90, cons_operator=">=")

    assert not at.exception
    assert any("70%" in str(i.value) for i in at.info)
    assert at.metric[0].value == ">= 70%"
    assert _run_button(at).disabled is False


def test_objective_contradicting_the_csvs_blocks_the_run(results_root):
    """'Únicos (<)' over CSVs computed at '>=' used to silently relabel the output."""
    at = _open_step5(results_root, threshold=70, cons_operator=">=")
    next(r for r in at.radio if r.label == "Objetivo").set_value("unique").run()

    assert any("não é o mesmo" in str(e.value) for e in at.error)
    assert _run_button(at).disabled is True


def _strip_threshold(results_root):
    """Replace the fixture's CSVs with one that declares no threshold anywhere."""
    csv_dir = results_root / "conservancy_csv"
    for stale in csv_dir.iterdir():
        stale.unlink()
    (csv_dir / "sp.csv").write_text(
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches,Minimum identity,Maximum identity,View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),80.00%,100.00%,d\n")
    return csv_dir


def test_unverifiable_csvs_warn_but_do_not_block(results_root):
    """Hand-made CSVs carrying no threshold anywhere stay usable, with the caveat stated."""
    _strip_threshold(results_root)

    at = _open_step5(results_root, threshold=70, cons_operator=">=")
    assert any("sem conseguir validá-lo" in str(w.value) for w in at.warning)
    assert _run_button(at).disabled is False


def test_unverifiable_csvs_let_the_user_set_the_threshold(results_root):
    """
    With nothing to read the threshold from, it must be asked — not defaulted.

    This branch used to take ``ss.threshold`` silently. That is the step-4 value only when step 4
    ran in the same session; after resuming an analysis from disk it was the initial 70, so POA2
    applied 70 to CSVs built at something else and said so in the metric.
    """
    _strip_threshold(results_root)
    at = _open_step5(results_root, threshold=70, cons_operator=">=")

    widget = next(n for n in at.number_input if "Limiar de identidade (t" in n.label)
    at = widget.set_value(90).run()

    assert not at.exception
    assert at.metric[0].value == ">= 90%"
    assert at.session_state["threshold"] == 90


def test_resuming_adopts_the_threshold_recorded_with_the_csvs(results_root, monkeypatch):
    """
    Resuming an analysis must not leave the threshold at its initial default.

    The cross-environment route (tune on Windows, run POA2 under WSL) reopens the analysis in a
    session that never ran step 4. ``_resume_from_results`` restored the reference FASTA and the
    POA1 table but not the threshold, so the session kept the 70 from ``setdefault`` while the
    CSVs on disk had been produced at another value.
    """
    cons_parser.write_metadata(str(results_root / "conservancy_csv"), 90, ">=", source="local")
    (results_root / "inputs").mkdir(exist_ok=True)
    (results_root / "inputs" / "proteins_all.fasta").write_text(">E_SP_p1\nGGAYIAMGG\n")

    at = AppTest.from_file(APP, default_timeout=90)
    at.session_state["run_path"] = str(results_root)
    at.session_state["conservancy_ready"] = False      # nothing done in this session
    at.run()
    at.sidebar.radio[0].set_value(STEP5).run()
    next(b for b in at.button if "Retomar" in b.label).click().run()

    assert not at.exception
    assert at.session_state["threshold"] == 90
    assert at.session_state["cons_operator"] == ">="
