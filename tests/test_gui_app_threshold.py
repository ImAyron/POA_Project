"""The threshold guard as the user meets it, driving the real Streamlit app headlessly.

This logic lives in ``poa/gui/app.py`` and is the user-facing half of the fix: POA2 must take the
threshold from the data rather than from the form, and must not let a run proceed under an
objective the CSVs were not computed with.
"""
import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from pathlib import Path  # noqa: E402

from streamlit.testing.v1 import AppTest  # noqa: E402

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
    epitopes.mkdir(parents=True)
    (epitopes / "SP_epitopes.fasta").write_text(">SP_E_BEPIPRED_1_5\nAYIAM\n")
    cc.run_conservancy_for_dir(
        str(epitopes),
        [("E_SP_p1", "GGAYIAMGG"), ("E_SP_p2", "TTAYIAMTT"), ("E_SP_p3", "GGAYIALGG")],
        threshold=70, out_dir=str(tmp_path / "conservancy_csv"), operator=">=")
    return tmp_path


def _open_step5(**state):
    at = AppTest.from_file(APP, default_timeout=90)
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
    at = _open_step5(threshold=90, cons_operator=">=")

    assert not at.exception
    assert any("70%" in str(i.value) for i in at.info)
    assert at.metric[0].value == ">= 70%"
    assert _run_button(at).disabled is False


def test_objective_contradicting_the_csvs_blocks_the_run(results_root):
    """'Únicos (<)' over CSVs computed at '>=' used to silently relabel the output."""
    at = _open_step5(threshold=70, cons_operator=">=")
    next(r for r in at.radio if r.label == "Objetivo").set_value("unique").run()

    assert any("não é o mesmo" in str(e.value) for e in at.error)
    assert _run_button(at).disabled is True


def test_unverifiable_csvs_warn_but_do_not_block(results_root):
    """Hand-made CSVs carrying no threshold anywhere stay usable, with the caveat stated."""
    csv_dir = results_root / "conservancy_csv"
    for stale in csv_dir.iterdir():
        stale.unlink()
    (csv_dir / "sp.csv").write_text(
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches,Minimum identity,Maximum identity,View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),80.00%,100.00%,d\n")

    at = _open_step5(threshold=70, cons_operator=">=")
    assert any("sem conseguir validá-lo" in str(w.value) for w in at.warning)
    assert _run_button(at).disabled is False
