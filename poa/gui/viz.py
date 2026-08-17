"""Visualization helpers for the GUI: a 2D epitope map (Plotly) and a 3D structure viewer.

Kept free of Streamlit imports so it can be unit-tested. The 3D viewer uses py3Dmol (3Dmol.js);
its generated HTML references 3Dmol.js from a CDN, so the 3D view needs internet at render time.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import plotly.express as px
import plotly.graph_objects as go


def _to_int(value, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (ValueError, TypeError):
        return default


def epitope_map_figure(df) -> go.Figure:
    """
    Build a 2D "feature map": each epitope drawn as a horizontal segment along the residue axis,
    one track per Protein_Specie, colored by prediction Method.

    Expects the standardized POA1 columns (Method, Specie, Protein, Initial Position,
    Final Position, Peptide Sequence).
    """
    fig = go.Figure()
    if df is None or df.empty:
        fig.update_layout(title="Sem epítopos para exibir")
        return fig

    d = df.copy()
    d["start"] = d["Initial Position"].map(_to_int)
    d["end"] = d["Final Position"].map(_to_int)
    d["length"] = (d["end"] - d["start"] + 1).clip(lower=1)
    d["track"] = d["Protein"].astype(str) + "_" + d["Specie"].astype(str)

    methods = list(d["Method"].unique())
    palette = px.colors.qualitative.Safe
    color_map = {m: palette[i % len(palette)] for i, m in enumerate(methods)}

    for m in methods:
        sub = d[d["Method"] == m]
        fig.add_bar(
            y=sub["track"],
            x=sub["length"],
            base=sub["start"],
            orientation="h",
            name=str(m),
            marker_color=color_map[m],
            opacity=0.7,
            customdata=sub[["Peptide Sequence", "start", "end"]].to_numpy(),
            hovertemplate="%{y}<br>%{customdata[0]}<br>pos %{customdata[1]}–%{customdata[2]}<extra>" + str(m) + "</extra>",
        )

    n_tracks = d["track"].nunique()
    fig.update_layout(
        barmode="overlay",
        xaxis_title="Posição do resíduo",
        yaxis_title="Proteína_Espécie",
        legend_title="Método",
        height=max(300, 40 * n_tracks + 140),
        margin=dict(l=10, r=10, t=40, b=10),
    )
    return fig


def epitope_ranges(df, protein: Optional[str] = None, specie: Optional[str] = None) -> List[Tuple[int, int]]:
    """Return [(start, end), ...] epitope residue ranges, optionally filtered by protein/specie."""
    if df is None or df.empty:
        return []
    d = df
    if protein is not None:
        d = d[d["Protein"].astype(str).str.upper() == protein.upper()]
    if specie is not None:
        d = d[d["Specie"].astype(str).str.upper() == specie.upper()]
    ranges = []
    for _, row in d.iterrows():
        start = _to_int(row["Initial Position"])
        end = _to_int(row["Final Position"])
        if start and end and end >= start:
            ranges.append((start, end))
    return ranges


def build_3dmol_view_html(
    pdb_text: str,
    ranges: List[Tuple[int, int]],
    chain: Optional[str] = None,
    base_style: str = "cartoon",
    show_surface: bool = False,
    epitope_color: str = "red",
    width: int = 760,
    height: int = 480,
) -> str:
    """
    Return self-rendering HTML (py3Dmol / 3Dmol.js) showing the structure with the given epitope
    residue ranges highlighted. Needs internet at render time (3Dmol.js is loaded from a CDN).
    """
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    view.setStyle({base_style: {"color": "lightgray"}})

    for start, end in ranges:
        selection = {"resi": list(range(start, end + 1))}
        if chain:
            selection["chain"] = chain
        view.addStyle(selection, {base_style: {"color": epitope_color}})

    if show_surface:
        view.addSurface(py3Dmol.VDW, {"opacity": 0.5, "color": "white"})

    view.zoomTo()
    return view._make_html()
