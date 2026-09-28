"""The analysis picker as the user meets it: create a run, reopen one, never mix two."""
import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from datetime import date  # noqa: E402
from pathlib import Path  # noqa: E402

from streamlit.testing.v1 import AppTest  # noqa: E402

from poa.gui import backend  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")


@pytest.fixture
def results_root(tmp_path, monkeypatch):
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    return tmp_path


def _start(**state):
    at = AppTest.from_file(APP, default_timeout=90)
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    return at


def test_no_analysis_selected_shows_the_picker_and_stops(results_root):
    at = _start()

    assert not at.exception
    assert any("Escolha a análise" in str(h.value) for h in at.header)
    assert not at.sidebar.radio            # the step selector is not reachable yet


def test_creating_an_analysis_binds_the_session_to_its_folder(results_root):
    at = _start()
    at.text_input[0].set_value("Teste DENV").run()
    next(b for b in at.button if "Criar" in b.label).click().run()

    assert not at.exception
    expected = results_root / backend.run_dir_name("Teste DENV")
    assert at.session_state["run_path"] == str(expected)
    assert expected.is_dir()
    assert at.sidebar.radio               # the steps are now reachable
    assert any(expected.name in str(c.value) for c in at.sidebar.code)


def test_existing_analyses_are_offered(results_root):
    backend.create_run(results_root, "Teste A", date(2026, 1, 5))
    backend.create_run(results_root, "Teste B", date(2026, 9, 27))

    at = _start()
    options = at.selectbox[0].options if at.selectbox else []
    assert len(options) == 2
    # newest first, labelled with the name the user typed and how far it got
    assert options[0].startswith("Teste B — 27/09/2026")
    assert options[1].startswith("Teste A — 05/01/2026")


def test_every_step_renders_inside_an_analysis(results_root):
    """ctx is now bound by the picker rather than at import; no step may break on that."""
    ctx = backend.create_run(results_root, "Teste A", date(2026, 9, 27))
    at = _start(run_path=str(ctx.root))
    steps = at.sidebar.radio[0].options

    assert len(steps) == 7
    for step in steps:
        at.sidebar.radio[0].set_value(step).run()
        assert not at.exception, f"{step}: {at.exception}"


def test_switching_analysis_clears_the_previous_one(results_root):
    """State points at files inside a folder; carrying it over would write into the wrong run."""
    first = backend.create_run(results_root, "Teste A", date(2026, 9, 27))
    second = backend.create_run(results_root, "Teste B", date(2026, 9, 27))

    at = _start(run_path=str(first.root), proteins_path=str(first.inputs_dir / "p.fasta"),
                conservancy_ready=True)
    assert at.session_state["proteins_path"]

    next(b for b in at.sidebar.button if "Trocar" in b.label).click().run()
    at.selectbox[0].set_value(1).run()           # index 1: the other analysis
    next(b for b in at.button if "Abrir" in b.label).click().run()

    assert at.session_state["run_path"] in (str(first.root), str(second.root))
    assert at.session_state["proteins_path"] is None
    assert at.session_state["conservancy_ready"] is False
