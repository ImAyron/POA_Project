"""Tests for the GUI visualization helpers (2D epitope map + 3D viewer HTML)."""
import pandas as pd
import plotly.graph_objects as go
import pytest

from poa.gui import viz


def _df():
    return pd.DataFrame({
        "Method": ["Bepipred3.0", "Bepipred3.0", "NetCTL"],
        "Specie": ["SARS", "SARS", "SARS"],
        "Protein": ["SPIKE", "SPIKE", "SPIKE"],
        "ID_Sequence": ["NP1", "NP1", "-"],
        "Initial Position": [4, 12, 20],
        "Final Position": [8, 15, 28],
        "Peptide Sequence": ["AYIAM", "LMNK", "IAMKGVLMN"],
    })


def test_epitope_map_figure_has_trace_per_method():
    fig = viz.epitope_map_figure(_df())
    assert isinstance(fig, go.Figure)
    # one bar trace per distinct method
    names = {t.name for t in fig.data}
    assert names == {"Bepipred3.0", "NetCTL"}


def test_epitope_map_figure_empty():
    fig = viz.epitope_map_figure(pd.DataFrame())
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 0


def test_epitope_ranges_parses_and_filters():
    df = _df()
    all_ranges = viz.epitope_ranges(df, protein="SPIKE", specie="SARS")
    assert (4, 8) in all_ranges and (12, 15) in all_ranges and (20, 28) in all_ranges
    # filtering by a different species yields nothing
    assert viz.epitope_ranges(df, specie="MERS") == []


def test_build_3dmol_view_html():
    pdb = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    html = viz.build_3dmol_view_html(pdb, [(1, 2)], base_style="cartoon")
    assert isinstance(html, str) and len(html) > 100
    assert "3dmol" in html.lower() or "viewer" in html.lower()


def test_missing_py3dmol_explains_how_to_install_it(monkeypatch):
    """
    The 3D viewer's dependency is imported lazily, so it fails at render time.

    A bare ImportError reached the user as "No module named 'py3Dmol'", which says nothing about
    what to do or that only step 7 is affected. The import stays deferred — the rest of the GUI,
    the CLI and these tests run without py3Dmol — so the message is where the fix has to be.
    """
    import builtins

    real_import = builtins.__import__

    def without_py3dmol(name, *args, **kwargs):
        if name == "py3Dmol":
            raise ImportError("No module named 'py3Dmol'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_py3dmol)

    with pytest.raises(RuntimeError, match="pip install py3Dmol"):
        viz.build_3dmol_view_html("ATOM\nEND\n", [(1, 2)])
