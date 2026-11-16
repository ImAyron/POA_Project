"""Visualization helpers for the GUI: a 2D epitope map (Plotly) and a 3D structure viewer.

Kept free of Streamlit imports so it can be unit-tested. The 3D viewer is a self-contained HTML
page driving 3Dmol.js directly; the page loads 3Dmol.js from a CDN, so the 3D view needs internet
at render time and nothing else — there is no Python package behind it.

It used to go through py3Dmol, which only translates Python calls into the same JavaScript. That
was enough to paint residues, and not enough for anything the page has to *do* afterwards: a
clipping control, a legend that zooms to the region clicked, a tooltip naming the residue under
the pointer all need code that runs in the browser and holds the viewer object, which py3Dmol
keeps in a variable of its own naming.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from poa.core import structure_map

#: The 3Dmol.js build the viewer page loads. Pinned to the release the page was written against:
#: it calls ``setSlab``, ``addLabel``, ``setHoverable`` and ``pngURI`` by name, and an unpinned
#: URL would hand a future major version to code that was never run against it.
THREEDMOL_JS = "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"

#: How the main chain can be drawn.
STYLES = ("cartoon", "stick", "sphere")

#: What ``crop`` accepts in :func:`build_3dmol_view_html`.
CROP_MODES = ("all", "chain", "range", "epitopes")

#: What ``color_by`` accepts in :func:`color_regions`.
COLOR_MODES = ("region", "method", "single")

BASE_COLOR = "#d3d3d3"            # the part of the structure that is not an epitope
MULTI_METHOD_COLOR = "#4d4d4d"    # a region predicted by more than one method, when coloring by method

#: One color per region, cycled. Chosen to stay apart from each other and from the gray base under
#: the common color-vision deficiencies (Okabe-Ito, minus its yellow — unreadable on white — and
#: extended with Tol's darker hues so ten neighbours differ before the cycle repeats).
REGION_PALETTE = ("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00",
                  "#56B4E9", "#882255", "#117733", "#332288", "#AA4499")

#: Above this many regions the 3D labels start off hidden: they are drawn in front of the
#: structure, so past a few dozen they cover the thing they are labelling.
MAX_AUTO_LABELS = 40


def _to_int(value, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (ValueError, TypeError):
        return default


def method_colors(methods: Iterable[str]) -> Dict[str, str]:
    """
    One color per prediction method, in the order the methods are given.

    Shared by the 2D map and the 3D viewer so that "colored by method" means the same color in
    both: reading the two views side by side is the reason to have both.
    """
    palette = px.colors.qualitative.Safe
    return {str(m): palette[i % len(palette)] for i, m in enumerate(dict.fromkeys(methods))}


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
    color_map = method_colors(methods)

    for m in methods:
        sub = d[d["Method"] == m]
        fig.add_bar(
            y=sub["track"],
            x=sub["length"],
            base=sub["start"],
            orientation="h",
            name=str(m),
            marker_color=color_map[str(m)],
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


# --------------------------------------------------------------------------- regions
@dataclass(frozen=True)
class Epitope:
    """One distinct epitope range of a protein, with every method that predicted it."""

    start: int
    end: int
    methods: Tuple[str, ...] = ()
    peptide: str = ""


def epitope_list(df, protein: Optional[str] = None, specie: Optional[str] = None,
                 methods: Optional[Sequence[str]] = None) -> List[Epitope]:
    """
    The distinct epitopes of one protein/species, in sequence order.

    :func:`epitope_ranges` answers "which residues", one tuple per table row. This answers "which
    epitopes": two methods predicting the same stretch are one entry carrying both names, because
    on a structure they are one stretch of residues, and listing it twice would label it twice.

    Args:
        df: the standardized POA1 table.
        protein, specie: keep only this track, compared case-insensitively.
        methods: keep only these prediction methods; all of them when omitted.

    Returns:
        One :class:`Epitope` per distinct ``(start, end)``, sorted by position.
    """
    if df is None or df.empty:
        return []
    d = df
    if protein is not None:
        d = d[d["Protein"].astype(str).str.upper() == protein.upper()]
    if specie is not None:
        d = d[d["Specie"].astype(str).str.upper() == specie.upper()]
    if methods is not None:
        d = d[d["Method"].astype(str).isin([str(m) for m in methods])]

    found: Dict[Tuple[int, int], Dict[str, object]] = {}
    for _, row in d.iterrows():
        start = _to_int(row["Initial Position"])
        end = _to_int(row["Final Position"])
        if not (start and end and end >= start):
            continue
        entry = found.setdefault((start, end), {"methods": [], "peptide": ""})
        method = str(row["Method"])
        if method not in entry["methods"]:
            entry["methods"].append(method)
        peptide = row.get("Peptide Sequence", "")
        if not entry["peptide"] and isinstance(peptide, str):    # a missing cell is NaN, not ""
            entry["peptide"] = peptide
    return [Epitope(start, end, tuple(entry["methods"]), str(entry["peptide"]))
            for (start, end), entry in sorted(found.items())]


@dataclass
class Region:
    """
    A stretch of the structure that gets highlighted, and what it stands for.

    Carries both numberings on purpose. ``seq_*`` is where the stretch is in the sequence the
    predictions ran on — the numbers in the POA1 table — and ``pdb_*`` is what the structure file
    calls the same residues, which is what the viewer has to be told. Showing only one of them is
    how a figure and a table end up describing "residues 482-495" and "residues 1-14" with nobody
    able to tell they are the same epitope.
    """

    name: str
    seq_start: int
    seq_end: int
    pdb_start: int
    pdb_end: int
    methods: Tuple[str, ...] = ()
    epitopes: int = 1
    peptide: str = ""
    partial: bool = False          # some of its residues are missing from the structure file
    color: str = "red"


def build_regions(epitopes: Sequence[Epitope],
                  placed: Mapping[Tuple[int, int], Tuple[int, int]], *,
                  merge: bool = True, partial: Iterable[Tuple[int, int]] = (),
                  sequence: str = "") -> List[Region]:
    """
    Turn mapped epitopes into the regions the viewer labels.

    A real run does not produce a handful of epitopes. The MHC-II predictions of one polyprotein
    in ``results/`` are 668 overlapping 15-mers; labelled one by one they are 668 tags over a
    structure that is almost entirely colored, and no tag can be matched to a stretch. Merged,
    the same predictions are the few dozen contiguous stretches a person actually sees on the
    structure — so that is the default, and the count each region stands for is kept.

    Args:
        epitopes: what :func:`epitope_list` returned.
        placed: ``(seq_start, seq_end) -> (pdb_start, pdb_end)`` for the epitopes the structure
            has, i.e. ``dict(zip(report.sources, report.ranges))`` of a mapping report. An epitope
            absent from it is not in the structure and produces no region.
        merge: fuse epitopes that overlap or touch into one region. When false every epitope is
            its own region.
        partial: sequence ranges whose residues are only partly present in the file.
        sequence: the reference sequence, used to spell out a merged region's residues.

    Returns:
        Regions in sequence order, named ``R1..`` when merged and ``E1..`` otherwise. Colors are
        not assigned here; see :func:`color_regions`.
    """
    gapped = set(partial)
    regions: List[Region] = []
    for epitope in sorted(epitopes, key=lambda e: (e.start, e.end)):
        key = (epitope.start, epitope.end)
        if key not in placed:
            continue
        pdb_start, pdb_end = placed[key]
        last = regions[-1] if regions else None
        if merge and last is not None and epitope.start <= last.seq_end + 1:
            if epitope.end > last.seq_end:
                last.seq_end, last.pdb_end = epitope.end, pdb_end
            last.methods = last.methods + tuple(m for m in epitope.methods if m not in last.methods)
            last.epitopes += 1
            last.partial = last.partial or key in gapped
            continue
        regions.append(Region(name="", seq_start=epitope.start, seq_end=epitope.end,
                              pdb_start=pdb_start, pdb_end=pdb_end, methods=epitope.methods,
                              peptide=epitope.peptide, partial=key in gapped))

    prefix = "R" if merge else "E"
    for index, region in enumerate(regions, start=1):
        region.name = f"{prefix}{index}"
        if region.epitopes > 1:
            # Several peptides, so no single one describes it; the stretch itself does.
            region.peptide = sequence[region.seq_start - 1:region.seq_end] if sequence else ""
    return regions


def color_regions(regions: Sequence[Region], color_by: str = "region",
                  colors: Optional[Mapping[str, str]] = None,
                  single_color: str = "red") -> Sequence[Region]:
    """
    Assign each region its color, in place.

    Args:
        regions: what :func:`build_regions` returned.
        color_by: ``'region'`` gives every region its own color, so a label, a legend row and a
            stretch of structure can be matched by eye; ``'method'`` uses the 2D map's colors
            (pass :func:`method_colors` as ``colors``), with one neutral color for a region that
            more than one method predicted; ``'single'`` paints them all ``single_color``.
        colors: method name -> color, for ``color_by='method'``.
        single_color: the color for ``color_by='single'``.

    Returns:
        The same sequence, for chaining.
    """
    if color_by not in COLOR_MODES:
        raise ValueError(f"color_by must be one of {COLOR_MODES}, not {color_by!r}")
    colors = colors or {}
    for index, region in enumerate(regions):
        if color_by == "single":
            region.color = single_color
        elif color_by == "method":
            region.color = (colors.get(region.methods[0], single_color)
                            if len(region.methods) == 1
                            else MULTI_METHOD_COLOR if region.methods else single_color)
        else:
            region.color = REGION_PALETTE[index % len(REGION_PALETTE)]
    return regions


def regions_frame(regions: Sequence[Region], chain: str = "") -> pd.DataFrame:
    """
    The highlighted regions as a table: the written counterpart of the labels on the structure.

    The names in the first column are the ones drawn on the structure and listed in the viewer's
    legend, so a region can be looked up in either direction.
    """
    prefix = f"{chain}:" if chain else ""
    rows = []
    for region in regions:
        rows.append({
            "Região": region.name,
            "Sequência": f"{region.seq_start}–{region.seq_end}",
            "PDB": f"{prefix}{region.pdb_start}–{region.pdb_end}",
            "Resíduos": region.seq_end - region.seq_start + 1,
            "Epítopos": region.epitopes,
            "Métodos": ", ".join(region.methods),
            "Peptídeo": region.peptide,
            "Observação": "parcialmente ausente no cristal" if region.partial else "",
        })
    return pd.DataFrame(rows, columns=["Região", "Sequência", "PDB", "Resíduos", "Epítopos",
                                       "Métodos", "Peptídeo", "Observação"])


# --------------------------------------------------------------------------- 3D viewer
_NAMED_COLORS = {"red": (255, 0, 0), "white": (255, 255, 255), "black": (0, 0, 0),
                 "lightgray": (211, 211, 211), "lightgrey": (211, 211, 211)}


def _rgb(color: str) -> Optional[Tuple[int, int, int]]:
    """``'#rrggbb'``, ``'rgb(r, g, b)'`` or one of a few names -> ``(r, g, b)``; None otherwise."""
    text = str(color).strip().lower()
    if re.fullmatch(r"#[0-9a-f]{6}", text):
        return tuple(int(text[i:i + 2], 16) for i in (1, 3, 5))
    match = re.fullmatch(r"rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", text)
    if match:
        return tuple(min(255, int(part)) for part in match.groups())
    return _NAMED_COLORS.get(text)


def _text_color(background: str) -> str:
    """Black or white, whichever reads better on ``background`` (white when it cannot be told)."""
    rgb = _rgb(background)
    if rgb is None:
        return "#ffffff"
    red, green, blue = rgb
    return "#000000" if (0.299 * red + 0.587 * green + 0.114 * blue) > 150 else "#ffffff"


def _residue_selection(start: int, end: int, chain) -> Dict[str, object]:
    """
    A 3Dmol.js selection for residues ``start..end``.

    Written as the range string 3Dmol.js understands (``'482-495'``) rather than the list of every
    number: it tests a list by walking it once per atom, and a merged region can span thousands of
    residues of a structure with tens of thousands of atoms. Negative numbers are fine in a range
    — the library's pattern is ``(-?\\d+)-(-?\\d+)`` — and they occur: chain A of PDB 4C23 starts
    at -5.

    ``chain`` is one identifier or a list of them (3Dmol.js matches any of a list); falsy means
    every chain.
    """
    selection: Dict[str, object] = {"resi": f"{start}-{end}"}
    if chain:
        selection["chain"] = chain
    return selection


def _visible_selection(crop: str, chain: Optional[str], crop_range: Optional[Tuple[int, int]],
                       crop_radius: float, region_selections: List[dict]) -> Optional[dict]:
    """The part of the structure to keep on screen, as a 3Dmol.js selection; None for all of it."""
    if crop == "chain" and chain:
        return {"chain": chain}
    if crop == "range" and crop_range:
        low, high = sorted(int(v) for v in crop_range)
        return _residue_selection(low, high, chain)
    if crop == "epitopes" and region_selections:
        union = {"or": region_selections}
        if crop_radius and crop_radius > 0:
            # Whole residues: an atom cut off from the rest of its side chain reads as damage.
            return {"within": {"distance": float(crop_radius), "sel": union}, "byres": True}
        return union
    return None


def _script_json(value) -> str:
    """
    JSON that is safe to place inside a ``<script>`` element.

    The structure file is somebody's upload and is embedded whole. ``json.dumps`` escapes quotes
    and newlines but leaves ``</script>`` intact, which would end the script element in the middle
    of a string and hand the rest of the file to the HTML parser — in a frame Streamlit runs with
    access to the app. Escaping every ``<`` (and ``>``, ``&``) as ``\\uXXXX`` keeps the text
    identical once parsed and leaves the HTML parser nothing to recognise.
    """
    return (json.dumps(value, ensure_ascii=True)
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026"))


def build_3dmol_view_html(
    pdb_text: str,
    ranges: Sequence[Tuple[int, int]],
    chain: Optional[str] = None,
    base_style: str = "cartoon",
    show_surface: bool = False,
    epitope_color: str = "red",
    height: int = 560,
    *,
    regions: Optional[Sequence[Region]] = None,
    copies: Sequence[str] = (),
    show_labels: Optional[bool] = None,
    side_chains: bool = False,
    crop: str = "all",
    crop_range: Optional[Tuple[int, int]] = None,
    crop_radius: float = 8.0,
    title: str = "",
) -> str:
    """
    Return a self-rendering HTML page (3Dmol.js) showing the structure with regions highlighted.

    The page is the whole viewer: the structure, a toolbar (fit, spin, labels, background, PNG and
    two sliders that clip the structure in depth), a tooltip naming the residue under the pointer,
    and a legend listing every highlighted region that zooms to the one clicked.

    One runtime requirement is not this function's to satisfy: **internet**, because the page
    loads 3Dmol.js from a CDN. The fetch happens in the viewer's browser, so nothing can be
    checked here; the page itself says so when the library does not arrive, instead of staying
    blank.

    Args:
        pdb_text: contents of the structure file, PDB or mmCIF (told apart by content).
        ranges: residue ranges to highlight, ``(start, end)`` inclusive, **in the file's own
            numbering** (see :mod:`poa.core.structure_map`).
        chain: restrict the highlight to this chain; every chain when omitted.
        base_style: one of :data:`STYLES`.
        show_surface: draw the molecular surface, colored like the residues under it.
        epitope_color: the highlight color when ``regions`` is not given.
        height: height of the whole page in pixels; the width follows the container.
        regions: what each entry of ``ranges`` stands for, index for index — its name, color and
            description. Without it the ranges are named ``R1..`` and painted ``epitope_color``.
        copies: other chains that are copies of ``chain`` under the same numbering (see
            :func:`poa.core.structure_map.chain_copies`); the regions are highlighted on them
            too. Ignored without ``chain``, where every chain is highlighted already.
        show_labels: start with the region names drawn on the structure. By default they are on
            up to :data:`MAX_AUTO_LABELS` regions; the toolbar toggles them either way.
        side_chains: also draw the side chains of the highlighted residues (cartoon only).
        crop: one of :data:`CROP_MODES` — ``'all'``; ``'chain'`` keeps only ``chain``;
            ``'range'`` keeps ``crop_range``; ``'epitopes'`` keeps the highlighted regions and
            everything within ``crop_radius``.
        crop_range: ``(start, end)`` in the file's numbering, for ``crop='range'``.
        crop_radius: ångströms around the regions to keep, for ``crop='epitopes'``; 0 keeps the
            regions alone.
        title: names the PNG the toolbar saves.

    Raises:
        ValueError: an unknown ``base_style`` or ``crop``, or ``regions`` that does not describe
            ``ranges``.
    """
    if base_style not in STYLES:
        raise ValueError(f"base_style must be one of {STYLES}, not {base_style!r}")
    if crop not in CROP_MODES:
        raise ValueError(f"crop must be one of {CROP_MODES}, not {crop!r}")

    ranges = [(int(start), int(end)) for start, end in ranges]
    if regions is None:
        regions = [Region(name=f"R{index}", seq_start=start, seq_end=end, pdb_start=start,
                          pdb_end=end, color=epitope_color)
                   for index, (start, end) in enumerate(ranges, start=1)]
    elif [(r.pdb_start, r.pdb_end) for r in regions] != ranges:
        raise ValueError("regions must describe ranges one for one, in the same order")
    if show_labels is None:
        show_labels = len(regions) <= MAX_AUTO_LABELS

    prefix = f"{chain}:" if chain else ""
    chains = [chain] + [c for c in copies if c != chain] if chain else []
    crowded = len(regions) > 12
    entries, selections = [], []
    for region in regions:
        selection = _residue_selection(region.pdb_start, region.pdb_end,
                                       chains if len(chains) > 1 else chain)
        selections.append(selection)
        sequence_span = f"{region.seq_start}–{region.seq_end}"
        pdb_span = f"{prefix}{region.pdb_start}–{region.pdb_end}"
        same = (region.seq_start, region.seq_end) == (region.pdb_start, region.pdb_end)
        detail = [", ".join(region.methods)] if region.methods else []
        if region.epitopes > 1:
            detail.append(f"{region.epitopes} epítopos")
        if region.partial:
            detail.append("parcial no cristal")
        entries.append({
            "name": region.name,
            # The span rides on the label while there is room for it; past a dozen regions the
            # name alone is drawn and the span is one click away in the legend.
            "label": region.name if crowded else f"{region.name} · {sequence_span}",
            "color": region.color,
            "fg": _text_color(region.color),
            "chains": chains,
            "start": region.pdb_start,
            "end": region.pdb_end,
            "sel": selection,
            "where": (f"resíduos {pdb_span}" if same
                      else f"seq. {sequence_span} → PDB {pdb_span}"),
            "detail": " · ".join(detail),
            "peptide": region.peptide,
        })

    visible = _visible_selection(crop, chain, crop_range, crop_radius, selections)
    note = ""
    if crop == "epitopes" and visible is None:
        note = ("Não há região destacada para recortar em volta; a estrutura é mostrada inteira.")
    elif crop != "all" and visible is None:
        note = "O recorte escolhido não se aplica a esta estrutura; ela é mostrada inteira."

    config = {
        "pdb": pdb_text,
        # 3Dmol.js names its mmCIF parser 'mcif'; the file's own format decides, not its name.
        "format": "mcif" if structure_map.structure_format(pdb_text) == "cif" else "pdb",
        "style": base_style,
        "baseColor": BASE_COLOR,
        "surface": bool(show_surface),
        "labels": bool(show_labels),
        "sidechains": bool(side_chains),
        "regions": entries,
        "visible": visible,
        "note": note,
        "png": (re.sub(r"[^A-Za-z0-9_.-]+", "_", title).strip("_") or "estrutura") + ".png",
    }
    return (_VIEWER_PAGE
            .replace("__HEIGHT__", str(int(height)))
            .replace("__THREEDMOL_JS__", THREEDMOL_JS)
            .replace("__CONFIG__", _script_json(config)))


# The page. Plain placeholders rather than str.format: the CSS and the JavaScript are full of
# braces. Everything that comes from data reaches the DOM through textContent, never innerHTML.
_VIEWER_PAGE = r"""<!doctype html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<style>
  * { box-sizing: border-box; }
  html, body { margin: 0; padding: 0; background: transparent; color: #1f2430;
               font: 13px/1.35 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
  .poa3d { display: flex; flex-direction: column; height: __HEIGHT__px; overflow: hidden;
           border: 1px solid #d0d4da; border-radius: 8px; background: #fff; }
  .bar { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 8px; padding: 6px 8px;
         background: #f4f5f7; border-bottom: 1px solid #d0d4da; }
  .bar button { font: inherit; color: inherit; background: #fff; border: 1px solid #c3c8d0;
                border-radius: 6px; padding: 3px 9px; cursor: pointer; }
  .bar button:hover { background: #eaf1fb; border-color: #7aa5dc; }
  .bar button[aria-pressed="true"] { background: #dbe8fa; border-color: #4b82c9; }
  .bar .cut { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 10px;
              margin-left: auto; }
  .bar label { display: flex; align-items: center; gap: 5px; color: #4a5261; white-space: nowrap; }
  .bar input[type=range] { width: 92px; margin: 0; }
  .main { display: flex; flex: 1; min-height: 0; }
  #stage { position: relative; flex: 1; min-width: 0; }
  #viewer { position: absolute; top: 0; right: 0; bottom: 0; left: 0; }
  #msg { position: absolute; left: 8px; bottom: 8px; max-width: calc(100% - 16px); z-index: 5;
         padding: 5px 9px; border-radius: 6px; background: #eef2f7; color: #3c4553;
         border: 1px solid #d0d4da; }
  #msg.bad { background: #fdecec; color: #8a1c1c; border-color: #f0b4b4; }
  #tip { position: absolute; display: none; z-index: 6; max-width: 280px; pointer-events: none;
         padding: 6px 9px; border-radius: 6px; background: rgba(28, 32, 40, .93); color: #fff;
         box-shadow: 0 2px 8px rgba(0, 0, 0, .25); }
  #tip b { display: block; }
  #tip .hit { display: flex; align-items: baseline; gap: 6px; margin-top: 3px; }
  #tip .sub { color: #c9cfd9; }
  .sw { flex: none; display: inline-block; width: 11px; height: 11px; border-radius: 3px;
        border: 1px solid rgba(0, 0, 0, .25); }
  #legend { flex: none; width: 250px; overflow-y: auto; background: #fafbfc;
            border-left: 1px solid #d0d4da; }
  #legend h4 { margin: 0; padding: 8px 10px 2px; font-size: 13px; }
  #legend p { margin: 0; padding: 0 10px 6px; color: #5b6472; font-size: 12px; }
  #legend .row { display: block; width: 100%; text-align: left; font: inherit; color: inherit;
                 background: none; border: 0; border-top: 1px solid #e6e9ee; padding: 6px 10px;
                 cursor: pointer; }
  #legend .row:hover { background: #eef3fa; }
  #legend .row.on { background: #dbe8fa; }
  #legend .head { display: flex; align-items: baseline; gap: 6px; }
  #legend .nm { font-weight: 600; }
  #legend .wh { color: #3c4553; }
  #legend .dt { color: #5b6472; font-size: 12px; margin-left: 17px; }
  #legend .pp { margin-left: 17px; color: #5b6472; white-space: nowrap; overflow: hidden;
                text-overflow: ellipsis; font: 11px/1.4 ui-monospace, Consolas, monospace; }
  #legend .rest { display: flex; align-items: center; gap: 6px; padding: 7px 10px;
                  border-top: 1px solid #e6e9ee; color: #5b6472; font-size: 12px; }
  @media (max-width: 640px) {
    .main { flex-direction: column; }
    #legend { width: auto; max-height: 38%; border-left: 0; border-top: 1px solid #d0d4da; }
  }
</style>
</head>
<body>
<div class="poa3d">
  <div class="bar">
    <button type="button" id="b-fit" title="Volta a enquadrar toda a parte visível">Enquadrar</button>
    <button type="button" id="b-spin" aria-pressed="false" title="Gira a estrutura">Girar</button>
    <button type="button" id="b-labels" aria-pressed="false"
            title="Mostra ou esconde o nome de cada região sobre a estrutura">Rótulos</button>
    <button type="button" id="b-bg" aria-pressed="false" title="Alterna o fundo">Fundo escuro</button>
    <button type="button" id="b-png" title="Salva a vista atual como imagem">Salvar PNG</button>
    <div class="cut" title="Corta a estrutura em profundidade, perpendicular à tela. Gire para escolher onde o corte passa.">
      <label>Corte frontal <input type="range" id="s-near" min="0" max="95" value="0"></label>
      <label>Corte traseiro <input type="range" id="s-far" min="0" max="95" value="0"></label>
      <button type="button" id="b-uncut" title="Desfaz os dois cortes">Sem corte</button>
    </div>
  </div>
  <div class="main">
    <div id="stage">
      <div id="viewer"></div>
      <div id="tip"></div>
      <div id="msg">Carregando o visualizador 3D…</div>
    </div>
    <div id="legend"></div>
  </div>
</div>
<script>var CFG = __CONFIG__;</script>
<script src="__THREEDMOL_JS__"></script>
<script>
(function () {
  "use strict";
  var cfg = CFG;
  var msg = document.getElementById("msg");

  function say(text, bad) {
    msg.textContent = text || "";
    msg.className = bad ? "bad" : "";
    msg.style.display = text ? "block" : "none";
  }

  function el(tag, cls, text) {
    var node = document.createElement(tag);
    if (cls) { node.className = cls; }
    if (text !== undefined) { node.textContent = text; }
    return node;
  }

  if (typeof $3Dmol === "undefined") {
    say("Não foi possível carregar o 3Dmol.js. O visualizador 3D precisa de internet no " +
        "momento de desenhar; o mapa 2D não depende dele.", true);
    return;
  }
  try {
    start();
  } catch (err) {
    say("Não foi possível desenhar a estrutura: " + (err && err.message ? err.message : err), true);
  }

  function start() {
    var stage = document.getElementById("stage");
    var viewer = $3Dmol.createViewer(document.getElementById("viewer"), {backgroundColor: "white"});
    viewer.addModel(cfg.pdb, cfg.format);

    var all = viewer.selectedAtoms({});
    if (!all.length) {
      say("O arquivo não contém átomos que o visualizador consiga ler.", true);
      return;
    }

    // ---- what stays on screen (the crop chosen in the form)
    var keep = null, note = cfg.note || "";
    if (cfg.visible) {
      var chosen = viewer.selectedAtoms(cfg.visible);
      if (chosen.length) {
        keep = {};
        chosen.forEach(function (a) { keep[a.index] = true; });
      } else {
        note = "O recorte pedido não contém nenhum átomo desta estrutura; ela é mostrada inteira.";
      }
    }
    function shown(a) { return keep === null || keep[a.index] === true; }
    function hidden(a) { return !shown(a); }
    var visSel = keep === null ? {} : {predicate: shown};
    var visible = keep === null ? all : all.filter(shown);

    // ---- which region(s) an atom belongs to
    function inRegion(r, a) {
      return (!r.chains.length || r.chains.indexOf(a.chain) !== -1) &&
             a.resi >= r.start && a.resi <= r.end;
    }
    function regionsAt(a) {
      return cfg.regions.filter(function (r) { return inRegion(r, a); });
    }
    function colorAt(a) {
      var hits = regionsAt(a);               // the last one wins, as it does when painting
      return hits.length ? hits[hits.length - 1].color : "#f2f2f2";
    }

    // ---- styles
    function styleFor(color, region) {
      var style = {};
      style[cfg.style] = {color: color};
      if (region && cfg.sidechains && cfg.style === "cartoon") {
        style.stick = {color: color, radius: 0.18};
      }
      return style;
    }
    viewer.setStyle({}, styleFor(cfg.baseColor, false));
    cfg.regions.forEach(function (r) { viewer.setStyle(r.sel, styleFor(r.color, true)); });
    if (keep !== null) { viewer.setStyle({predicate: hidden}, {}); }

    if (cfg.surface) {
      var protein = {hetflag: false, predicate: shown};
      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.85, colorfunc: colorAt},
                        protein, protein);
    }

    // ---- labels drawn on the structure
    function anchor(atoms) {
      var x = 0, y = 0, z = 0;
      atoms.forEach(function (a) { x += a.x; y += a.y; z += a.z; });
      x /= atoms.length; y /= atoms.length; z /= atoms.length;
      // The centroid of a curved stretch floats off the chain; sit on the nearest atom instead.
      var best = atoms[0], bestD = Infinity;
      atoms.forEach(function (a) {
        var d = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y) + (a.z - z) * (a.z - z);
        if (d < bestD) { bestD = d; best = a; }
      });
      return {x: best.x, y: best.y, z: best.z};
    }
    var regionAtoms = cfg.regions.map(function (r) {
      return viewer.selectedAtoms(r.sel).filter(shown);
    });
    var labelsOn = !!cfg.labels;
    var labelsBtn = document.getElementById("b-labels");
    function drawLabels() {
      viewer.removeAllLabels();
      labelsBtn.setAttribute("aria-pressed", labelsOn ? "true" : "false");
      if (!labelsOn) { return; }
      cfg.regions.forEach(function (r, i) {
        // One label per chain the region is drawn on: a single label for a dimer would sit on
        // one copy, or between the two, and name neither.
        var perChain = {};
        regionAtoms[i].forEach(function (a) { (perChain[a.chain] = perChain[a.chain] || []).push(a); });
        Object.keys(perChain).forEach(function (c) {
          viewer.addLabel(r.label, {
            position: anchor(perChain[c]), backgroundColor: r.color, backgroundOpacity: 0.9,
            fontColor: r.fg, fontSize: 12, borderThickness: 0, inFront: true,
            showBackground: true, alignment: "center"
          });
        });
      });
    }

    // ---- tooltip: the residue under the pointer and the regions it is part of
    var tip = document.getElementById("tip");
    var mouse = {x: 0, y: 0};
    function placeTip() {
      var left = mouse.x + 14, top = mouse.y + 14;
      if (left + tip.offsetWidth > stage.clientWidth) { left = mouse.x - tip.offsetWidth - 10; }
      if (top + tip.offsetHeight > stage.clientHeight) { top = mouse.y - tip.offsetHeight - 10; }
      tip.style.left = Math.max(0, left) + "px";
      tip.style.top = Math.max(0, top) + "px";
    }
    stage.addEventListener("mousemove", function (e) {
      var box = stage.getBoundingClientRect();
      mouse.x = e.clientX - box.left;
      mouse.y = e.clientY - box.top;
      if (tip.style.display === "block") { placeTip(); }
    });
    stage.addEventListener("mouseleave", function () { tip.style.display = "none"; });
    function showTip(atom) {
      while (tip.firstChild) { tip.removeChild(tip.firstChild); }
      tip.appendChild(el("b", "", atom.resn + " " + atom.resi + (atom.icode || "").trim() +
                                  " · cadeia " + atom.chain));
      var hits = regionsAt(atom);
      if (!hits.length) { tip.appendChild(el("span", "sub", "fora das regiões destacadas")); }
      hits.forEach(function (r) {
        var line = el("div", "hit");
        var sw = el("span", "sw");
        sw.style.background = r.color;
        line.appendChild(sw);
        var text = el("span", "", r.name + " · " + r.where);
        if (r.detail) { text.appendChild(el("div", "sub", r.detail)); }
        line.appendChild(text);
        tip.appendChild(line);
      });
      tip.style.display = "block";
      placeTip();
    }
    // Every visible atom answers the pointer, unless there are so many that building a pick
    // target for each would stall the page; then only the highlighted ones do.
    var hoverSel = visible.length <= 40000 ? visSel : {predicate: function (a) {
      return shown(a) && regionsAt(a).length > 0;
    }};
    viewer.setHoverDuration(120);
    viewer.setHoverable(hoverSel, true, showTip, function () { tip.style.display = "none"; });

    // ---- framing and the depth cut
    var lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    visible.forEach(function (a) {
      lo[0] = Math.min(lo[0], a.x); lo[1] = Math.min(lo[1], a.y); lo[2] = Math.min(lo[2], a.z);
      hi[0] = Math.max(hi[0], a.x); hi[1] = Math.max(hi[1], a.y); hi[2] = Math.max(hi[2], a.z);
    });
    var depth = Math.max(10, Math.sqrt(Math.pow(hi[0] - lo[0], 2) + Math.pow(hi[1] - lo[1], 2) +
                                       Math.pow(hi[2] - lo[2], 2)));
    var nearS = document.getElementById("s-near"), farS = document.getElementById("s-far");
    var base = {near: -depth / 1.9, far: depth / 2};
    function applySlab() {
      var span = base.far - base.near;
      var near = base.near + span * Number(nearS.value) / 100;
      var far = base.far - span * Number(farS.value) / 100;
      if (far < near + 1) { far = near + 1; }
      viewer.setSlab(near, far);
      viewer.render();
    }
    // zoomTo() resets the slab to fit what it zooms to, which on a single region clips away the
    // rest of the structure. The slab is put back to span the whole visible part: around its own
    // centre for the full view, twice as deep around a region, whose centre can sit at the edge.
    function frame(sel, wide, ms) {
      viewer.zoomTo(sel, ms || 0);
      var k = wide ? 2 : 1;
      base = {near: -k * depth / 1.9, far: k * depth / 2};
      nearS.value = 0;
      farS.value = 0;
      applySlab();
    }
    nearS.addEventListener("input", applySlab);
    farS.addEventListener("input", applySlab);
    document.getElementById("b-uncut").addEventListener("click", function () {
      nearS.value = 0;
      farS.value = 0;
      applySlab();
    });

    // ---- legend: one row per region, click to zoom
    var legend = document.getElementById("legend");
    var active = null;
    function setActive(row) {
      if (active) { active.classList.remove("on"); }
      active = row;
      if (active) { active.classList.add("on"); }
    }
    legend.appendChild(el("h4", "", "Regiões destacadas (" + cfg.regions.length + ")"));
    legend.appendChild(el("p", "", cfg.regions.length
      ? "Clique em uma região para aproximar; clique de novo para voltar."
      : "Nenhuma região destacada nesta estrutura."));
    cfg.regions.forEach(function (r, i) {
      var row = el("button", "row");
      row.type = "button";
      var head = el("div", "head");
      var sw = el("span", "sw");
      sw.style.background = r.color;
      head.appendChild(sw);
      head.appendChild(el("span", "nm", r.name));
      head.appendChild(el("span", "wh", r.where));
      row.appendChild(head);
      if (r.detail) { row.appendChild(el("div", "dt", r.detail)); }
      if (r.peptide) {
        var pep = el("div", "pp", r.peptide);
        pep.title = r.peptide;
        row.appendChild(pep);
      }
      if (!regionAtoms[i].length) {
        row.appendChild(el("div", "dt", "fora da parte visível"));
      }
      row.addEventListener("click", function () {
        if (active === row || !regionAtoms[i].length) {
          setActive(null);
          frame(visSel, false, 400);
          return;
        }
        setActive(row);
        frame({predicate: function (a) { return shown(a) && inRegion(r, a); }}, true, 400);
      });
      legend.appendChild(row);
    });
    var rest = el("div", "rest");
    var gray = el("span", "sw");
    gray.style.background = cfg.baseColor;
    rest.appendChild(gray);
    rest.appendChild(el("span", "", "restante da estrutura"));
    legend.appendChild(rest);

    // ---- toolbar
    document.getElementById("b-fit").addEventListener("click", function () {
      setActive(null);
      frame(visSel, false, 400);
    });
    var spinning = false;
    document.getElementById("b-spin").addEventListener("click", function () {
      spinning = !spinning;
      this.setAttribute("aria-pressed", spinning ? "true" : "false");
      viewer.spin(spinning ? "y" : false);
    });
    labelsBtn.addEventListener("click", function () {
      labelsOn = !labelsOn;
      drawLabels();
      viewer.render();
    });
    var dark = false;
    document.getElementById("b-bg").addEventListener("click", function () {
      dark = !dark;
      this.setAttribute("aria-pressed", dark ? "true" : "false");
      viewer.setBackgroundColor(dark ? "black" : "white");
      viewer.render();
    });
    document.getElementById("b-png").addEventListener("click", function () {
      viewer.render();
      var link = document.createElement("a");
      link.href = viewer.pngURI();
      link.download = cfg.png;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    });
    window.addEventListener("resize", function () { viewer.resize(); });

    drawLabels();
    frame(visSel, false, 0);
    say(note, false);
    window.poaViewerReady = true;
  }
})();
</script>
</body>
</html>
"""


