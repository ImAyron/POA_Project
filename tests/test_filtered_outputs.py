"""Telling the POA2 selection apart from the conservancy input it was filtered out of.

The CSVs of step 4 hold *every* epitope at any threshold — the identity threshold decides what
the "percent of matches" column counts, it never drops a row (see
``conservancy_client.conservancy_rows``). The filtered selection is the POA2 spreadsheet. With
the only CSVs in an analysis being the unfiltered ones, tightening a filter and opening "the CSV"
showed the whole set, reading exactly like a filter that had been ignored. So step 4 says what
its files are, and step 6 offers the selection as a CSV too.
"""
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from streamlit.testing.v1 import AppTest  # noqa: E402

from poa.core import pipeline  # noqa: E402
from poa.services import conservancy_client as cc  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")
STEP4 = "4 · Conservancy Analysis"
STEP6 = "6 · Resultados"


@pytest.fixture
def results_root(tmp_path, monkeypatch):
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    epitopes = tmp_path / "poa1_out" / "Conservancy Analysis"
    epitopes.mkdir(parents=True)
    (epitopes / "SP_epitopes.fasta").write_text(">SP_E_BEPIPRED_1_5\nAYIAM\n")
    cc.run_conservancy_for_dir(
        str(epitopes),
        [("E_SP_p1", "GGAYIAMGG"), ("E_SP_p2", "TTAYIAMTT"), ("E_SP_p3", "GGAYIALGG")],
        threshold=70, out_dir=str(tmp_path / "conservancy_csv"), operator=">=")
    return tmp_path


def _open(root, step, **state):
    at = AppTest.from_file(APP, default_timeout=90)
    at.session_state["run_path"] = str(root)
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    return at.sidebar.radio[0].set_value(step).run()


# --------------------------------------------------------------------------- step 4
def test_step4_says_its_csvs_are_the_unfiltered_input(results_root):
    at = _open(results_root, STEP4, poa1_result=object())

    assert not at.exception
    said = " ".join(str(i.value) for i in at.info)
    assert "entrada do POA2" in said
    assert "não descarta linhas" in said


def test_step4_download_is_labelled_as_the_whole_set(results_root):
    at = _open(results_root, STEP4, poa1_result=object())

    assert any("todos os epítopos" in str(b.label) for b in at.download_button)


# --------------------------------------------------------------------------- step 6
@pytest.fixture
def poa2_result(results_root):
    """A POA2 selection on disk, as a finished run leaves it."""
    selection = pd.DataFrame({
        "Epitope name": ["SP_E_BEPIPRED_1_5"],
        "Epitope sequence": ["AYIAM"],
        "Epitope length": [5],
        "Minimum identity(%)": [66.67],
        "Maximum identity(%)": [100.0],
    })
    out = results_root / "poa2_out"
    out.mkdir(exist_ok=True)
    xlsx = out / "POA2_analysis_ge70.xlsx"
    selection.to_excel(xlsx, index=False)
    return pipeline.Poa2Result(selection, str(xlsx), None, ">=")


def test_step6_offers_the_selection_as_csv(results_root, poa2_result):
    at = _open(results_root, STEP6, poa2_result=poa2_result, conservancy_ready=True)

    assert not at.exception
    assert any("Epítopos selecionados (.csv)" in str(b.label) for b in at.download_button)


def test_the_csv_holds_exactly_the_selection(poa2_result):
    """The bytes the button hands over — the proto carries neither them nor the file name."""
    from io import BytesIO

    from poa.gui import backend

    written = pd.read_csv(BytesIO(backend.results_csv(poa2_result.results)),
                          encoding="utf-8-sig")

    pd.testing.assert_frame_equal(written, poa2_result.results)


def test_the_csv_opens_as_utf8_in_excel():
    """Accented names survive the round trip Excel on Windows puts a CSV through."""
    from io import BytesIO

    from poa.gui import backend

    frame = pd.DataFrame({"Epitope name": ["PROTEÍNA_ESPÉCIE_1_5"]})
    data = backend.results_csv(frame)

    assert data.startswith(b"\xef\xbb\xbf")     # the BOM Excel needs to pick UTF-8
    back = pd.read_csv(BytesIO(data), encoding="utf-8-sig")
    assert back["Epitope name"].iloc[0] == "PROTEÍNA_ESPÉCIE_1_5"


def test_the_csv_of_an_empty_selection_is_still_readable():
    from poa.gui import backend

    assert backend.results_csv(pd.DataFrame()).startswith(b"\xef\xbb\xbf")


def test_step6_separates_the_selection_from_the_step4_input(results_root, poa2_result):
    at = _open(results_root, STEP6, poa2_result=poa2_result, conservancy_ready=True)

    said = " ".join(str(c.value) for c in at.caption)
    assert "etapa 4 são a entrada" in said
