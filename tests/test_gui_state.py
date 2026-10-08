"""The state policy of the GUI: what a step shows must survive navigating away from it.

Streamlit reruns the script on every interaction and renders only the step selected in the
sidebar. A widget that is not rendered loses its state, so any input that is not mirrored into
session state is reset to its own default the next time its step is opened. These tests pin the
two halves of the policy stated at the top of ``poa/gui/app.py``: the values round-trip, and
every session key a step reads is declared before any step runs.
"""
import ast
from pathlib import Path

import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from streamlit.testing.v1 import AppTest  # noqa: E402

from poa.services import conservancy_client as cc  # noqa: E402
from poa.services import tmhmm_client  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")
STEPS_FIRST = "1 · Proteínas & Parâmetros"
STEP4 = "4 · Conservancy Analysis"
STEP5 = "5 · POA2 (conservação + topologia)"


@pytest.fixture
def results_root(tmp_path, monkeypatch):
    """A results tree whose conservancy CSVs let step 5 render its filter form."""
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    monkeypatch.setattr(tmhmm_client, "is_available", lambda: True)

    epitopes = tmp_path / "poa1_out" / "Conservancy Analysis"
    epitopes.mkdir(parents=True)
    (epitopes / "SP_epitopes.fasta").write_text(">SP_E_BEPIPRED_1_5\nAYIAM\n")
    cc.run_conservancy_for_dir(
        str(epitopes),
        [("E_SP_p1", "GGAYIAMGG"), ("E_SP_p2", "TTAYIAMTT"), ("E_SP_p3", "GGAYIALGG")],
        threshold=70, out_dir=str(tmp_path / "conservancy_csv"), operator=">=")
    return tmp_path


def _open_step5(root):
    at = AppTest.from_file(APP, default_timeout=90)
    at.session_state["run_path"] = str(root)      # skip the analysis picker
    at.session_state["conservancy_ready"] = True
    at.session_state["proteins_path"] = "proteins.fasta"
    at.run()
    return at.sidebar.radio[0].set_value(STEP5).run()


def _number(at, label):
    return next(n for n in at.number_input if label in n.label)


def _checkbox(at, label):
    return next(c for c in at.checkbox if label in c.label)


def test_poa2_filters_survive_navigating_to_another_step(results_root):
    """
    Leaving step 5 and coming back must not reset the filters.

    The POA2 filter widgets carried no session key, so opening step 4 to check a CSV silently
    reverted them to 60/100/60 — and these are the parameters that decide which epitopes are
    selected. The reset was invisible: the form simply showed its defaults again.
    """
    at = _open_step5(results_root)

    at = _number(at, "Identidade mín").set_value(80).run()
    at = _number(at, "match (-m)").set_value(95).run()
    at = _checkbox(at, "também como filtro").set_value(True).run()
    assert at.session_state["params_poa2"]["imin"] == 80

    at = at.sidebar.radio[0].set_value(STEP4).run()
    at = at.sidebar.radio[0].set_value(STEP5).run()

    assert not at.exception
    assert _number(at, "Identidade mín").value == 80
    assert _number(at, "match (-m)").value == 95
    assert _checkbox(at, "também como filtro").value is True
    assert at.session_state["params_poa2"]["m"] == 95


def test_poa2_form_is_mirrored_into_session_state(results_root):
    """The key of record must hold the whole form, not a subset of it."""
    at = _open_step5(results_root)

    saved = at.session_state["params_poa2"]
    assert set(saved) == {"objective", "t", "imin", "imax", "m", "rf", "idf"}


def test_objective_still_follows_the_csvs_not_the_saved_form(results_root):
    """
    The one field the policy must NOT restore.

    The objective has to match the operator the CSVs were produced with — POA2 blocks the run
    otherwise. Restoring it from the saved form would carry an objective chosen over one analysis
    onto CSVs built the other way round, so it is always derived from the data.
    """
    at = _open_step5(results_root)                       # fixture CSVs are '>='
    at = next(r for r in at.radio if r.label == "Objetivo").set_value("unique").run()
    assert any("não é o mesmo" in str(e.value) for e in at.error)

    at = at.sidebar.radio[0].set_value(STEP4).run()
    at = at.sidebar.radio[0].set_value(STEP5).run()

    assert not at.exception
    assert next(r for r in at.radio if r.label == "Objetivo").value == "conserved"


def test_new_analysis_button_clears_everything_and_unbinds_the_folder(results_root):
    """
    "Nova análise" is the one path that wipes the forms, and it must land on the picker.

    Switching analyses keeps the filter parameters on purpose. Starting over is the opposite
    intent, so this clears the form keys and the step selector too and leaves ``run_path`` unset,
    which is what makes a *new* folder be created — reusing the current one would leave its
    poa1_out/ and conservancy_csv/ on disk to mix into the next test.
    """
    at = _open_step5(results_root)
    at = _number(at, "Identidade mín").set_value(80).run()
    assert at.session_state["params_poa2"]["imin"] == 80
    assert at.session_state["step_radio"] == STEP5

    at = next(b for b in at.sidebar.button if "Nova análise" in b.label).click().run()

    assert not at.exception
    assert at.session_state["run_path"] is None          # back to the analysis picker
    assert at.session_state["params_poa2"] == {}         # the form was cleared, not carried over
    assert at.session_state.get("step_radio") in (None, STEPS_FIRST)
    assert any("Escolha a análise" in str(h.value) for h in at.header)


def test_switching_analyses_keeps_the_filter_parameters(results_root):
    """The narrower reset: another folder, same parameters — the documented difference."""
    at = _open_step5(results_root)
    at = _number(at, "Identidade mín").set_value(80).run()

    at = next(b for b in at.sidebar.button if "Trocar de análise" in b.label).click().run()

    assert not at.exception
    assert at.session_state["run_path"] is None
    assert at.session_state["params_poa2"]["imin"] == 80


def test_no_upload_is_read_without_being_saved_first(results_root):
    """
    Every uploaded file must reach disk, because the uploader itself does not survive navigation.

    ``st.file_uploader`` returns ``None`` on any rerun where its step is not rendered, so a file
    read straight from the widget is gone the moment the user opens another step. The step-7 PDB
    was the one upload doing that: the 3D view vanished and the structure had to be sent again.
    This reads the module rather than driving the UI, so it also covers uploaders added later.
    """
    tree = ast.parse(Path(APP).read_text(encoding="utf-8"))
    offenders = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "getvalue"):
            continue
        # save_upload(...) / import_bepipred2(...) take the bytes and write them out; a getvalue()
        # anywhere else is a file that only ever existed inside that one rerun.
        parent = next((p for p in ast.walk(tree)
                       if isinstance(p, ast.Call) and node in ast.walk(p) and p is not node), None)
        if parent is None or getattr(parent.func, "attr", "") not in {
                "save_upload", "import_bepipred2"}:
            offenders.append(ast.unparse(node))

    assert not offenders, f"uploads lidos sem serem gravados: {offenders}"


def test_every_session_key_is_declared_before_use():
    """
    The policy is only safe if every ``ss.<key>`` exists before a step reads it.

    A step reached straight from the sidebar runs without the steps before it, so a key created
    on step 1 and read on step 5 would raise ``AttributeError`` for that user. ``_init_state`` is
    the single place that declares them, and this fails when a new key skips it.
    """
    tree = ast.parse(Path(APP).read_text(encoding="utf-8"))
    used = {node.attr for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name) and node.value.id == "ss"}
    declared = {arg.value for node in ast.walk(tree)
                if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "setdefault"
                for arg in node.args[:1] if isinstance(arg, ast.Constant)}
    mapping_methods = {"setdefault", "get", "keys"}   # ss is a mapping, not only an attr bag

    assert not (used - declared - mapping_methods)
