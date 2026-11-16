"""Tests for the GUI visualization helpers (2D epitope map + 3D viewer HTML)."""
import json

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


PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
    "END\n"
)


def _config(html):
    """The data the viewer page was built with, read back out of the page."""
    start = html.index("var CFG = ") + len("var CFG = ")
    return json.loads(html[start:html.index(";</script>", start)])


def test_build_3dmol_view_html():
    html = viz.build_3dmol_view_html(PDB, [(1, 2)], base_style="cartoon")

    assert viz.THREEDMOL_JS in html                      # the library is loaded from the pinned URL
    config = _config(html)
    assert config["pdb"] == PDB
    assert config["style"] == "cartoon"
    assert [r["name"] for r in config["regions"]] == ["R1"]
    assert config["regions"][0]["sel"] == {"resi": "1-2"}
    assert config["visible"] is None                     # nothing cropped unless asked


def test_the_viewer_needs_no_python_package(monkeypatch):
    """
    The 3D viewer is a page, not a wrapper: building it must not import anything optional.

    It used to go through py3Dmol, imported lazily, so a missing package only showed at render
    time — as "No module named 'py3Dmol'" until that was translated into install instructions.
    The page now talks to 3Dmol.js directly, and the way to keep that failure from coming back
    is to keep the import from coming back.
    """
    import builtins

    real_import = builtins.__import__

    def without_py3dmol(name, *args, **kwargs):
        if name == "py3Dmol":
            raise ImportError("No module named 'py3Dmol'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_py3dmol)

    assert "$3Dmol.createViewer" in viz.build_3dmol_view_html(PDB, [(1, 2)])


def test_an_uploaded_file_cannot_break_out_of_the_page():
    """
    The structure is somebody's upload, embedded whole in a ``<script>`` element.

    Streamlit runs the page in a frame with access to the app, so a file that could end the
    script element early would be running its own markup there. ``json.dumps`` does not escape
    ``</script>``; the page must come out with exactly the script elements the template has, and
    the text must still read back unchanged.
    """
    hostile = PDB + "REMARK </script><script>alert(1)</script><!-- & -->\n"
    clean = viz.build_3dmol_view_html(PDB, [(1, 2)])
    html = viz.build_3dmol_view_html(hostile, [(1, 2)])

    assert html.count("</script>") == clean.count("</script>")
    assert "alert(1)</script>" not in html
    assert _config(html)["pdb"] == hostile


# --------------------------------------------------------------------------- regions
def test_epitope_list_is_one_entry_per_stretch():
    """Two methods predicting the same residues are one stretch of structure, named once."""
    df = pd.concat([_df(), pd.DataFrame({
        "Method": ["NetCTL"], "Specie": ["SARS"], "Protein": ["SPIKE"], "ID_Sequence": ["-"],
        "Initial Position": [4], "Final Position": [8], "Peptide Sequence": ["AYIAM"]})])

    epitopes = viz.epitope_list(df, protein="SPIKE", specie="SARS")

    assert [(e.start, e.end) for e in epitopes] == [(4, 8), (12, 15), (20, 28)]
    assert epitopes[0].methods == ("Bepipred3.0", "NetCTL")
    assert epitopes[0].peptide == "AYIAM"
    # the method filter drops the stretches only the other method predicted
    only = viz.epitope_list(df, protein="SPIKE", specie="SARS", methods=["NetCTL"])
    assert [(e.start, e.end) for e in only] == [(4, 8), (20, 28)]
    assert viz.epitope_list(df, specie="MERS") == []


def test_overlapping_epitopes_become_one_named_region():
    """
    Sliding-window predictors emit hundreds of near-identical peptides.

    One analysis in ``results/`` has 668 MHC-II 15-mers on a single polyprotein. A label per
    peptide is 668 labels on a structure that is nearly all colored; merged, they are the
    contiguous stretches a person can actually point at, each knowing how many it stands for.
    """
    epitopes = [viz.Epitope(1, 15, ("MHCII",), "A" * 15), viz.Epitope(5, 19, ("MHCII",), "B" * 15),
                viz.Epitope(20, 24, ("PAP",), "CCCCC"),          # touches the previous one
                viz.Epitope(40, 50, ("PAP",), "D" * 11)]
    placed = {(e.start, e.end): (e.start + 481, e.end + 481) for e in epitopes}
    sequence = "".join(chr(65 + i % 26) for i in range(60))

    regions = viz.build_regions(epitopes, placed, merge=True, partial=[(5, 19)],
                                sequence=sequence)

    assert [r.name for r in regions] == ["R1", "R2"]
    first, second = regions
    assert (first.seq_start, first.seq_end) == (1, 24)
    assert (first.pdb_start, first.pdb_end) == (482, 505)      # the file's numbering, not ours
    assert first.epitopes == 3 and first.methods == ("MHCII", "PAP")
    assert first.partial is True                               # one member is gapped
    assert first.peptide == sequence[0:24]                     # the stretch, not one member
    assert second.epitopes == 1 and second.peptide == "D" * 11


def test_ungrouped_regions_are_the_epitopes_themselves():
    epitopes = [viz.Epitope(1, 15, ("MHCII",), "A" * 15), viz.Epitope(5, 19, ("MHCII",), "B" * 15)]
    placed = {(1, 15): (1, 15), (5, 19): (5, 19)}

    regions = viz.build_regions(epitopes, placed, merge=False)

    assert [(r.name, r.seq_start, r.seq_end) for r in regions] == [("E1", 1, 15), ("E2", 5, 19)]


def test_an_epitope_the_structure_lacks_makes_no_region():
    """Only what the mapping placed is highlighted; the rest is reported, never guessed."""
    epitopes = [viz.Epitope(1, 5, ("M",), "AAAAA"), viz.Epitope(100, 105, ("M",), "BBBBBB")]

    regions = viz.build_regions(epitopes, {(1, 5): (10, 14)})

    assert [(r.seq_start, r.pdb_start) for r in regions] == [(1, 10)]


def test_color_modes():
    def fresh():
        return [viz.Region("R1", 1, 5, 1, 5, methods=("A",)),
                viz.Region("R2", 9, 12, 9, 12, methods=("B",)),
                viz.Region("R3", 20, 30, 20, 30, methods=("A", "B"))]

    by_region = viz.color_regions(fresh(), "region")
    assert len({r.color for r in by_region}) == 3              # told apart by eye

    palette = viz.method_colors(["A", "B"])
    by_method = viz.color_regions(fresh(), "method", palette)
    assert [r.color for r in by_method] == [palette["A"], palette["B"], viz.MULTI_METHOD_COLOR]

    assert {r.color for r in viz.color_regions(fresh(), "single")} == {"red"}
    with pytest.raises(ValueError, match="color_by"):
        viz.color_regions(fresh(), "rainbow")


def test_the_2d_map_and_the_3d_viewer_agree_on_method_colors():
    """Reading the two views side by side only works if a method is one color in both."""
    fig = viz.epitope_map_figure(_df())
    colors = viz.method_colors(_df()["Method"].unique())

    assert {t.name: t.marker.color for t in fig.data} == colors


def test_regions_frame_names_match_the_viewer():
    regions = viz.build_regions([viz.Epitope(1, 5, ("M",), "AAAAA")], {(1, 5): (482, 486)})
    viz.color_regions(regions)

    frame = viz.regions_frame(regions, chain="A")
    legend = _config(viz.build_3dmol_view_html(PDB, [(482, 486)], chain="A", regions=regions))

    assert list(frame["Região"]) == [r["name"] for r in legend["regions"]] == ["R1"]
    assert frame.loc[0, "Sequência"] == "1–5" and frame.loc[0, "PDB"] == "A:482–486"
    # the legend says both numberings too, since they differ here
    assert "1–5" in legend["regions"][0]["where"] and "A:482–486" in legend["regions"][0]["where"]


# --------------------------------------------------------------------------- viewer options
def test_regions_must_describe_the_ranges():
    regions = [viz.Region("R1", 1, 5, 1, 5)]

    with pytest.raises(ValueError, match="one for one"):
        viz.build_3dmol_view_html(PDB, [(1, 2)], regions=regions)
    with pytest.raises(ValueError, match="base_style"):
        viz.build_3dmol_view_html(PDB, [(1, 2)], base_style="ribbon")
    with pytest.raises(ValueError, match="crop"):
        viz.build_3dmol_view_html(PDB, [(1, 2)], crop="half")


def test_negative_residue_numbers_survive_the_selection():
    """Expression tags push numbering below 1 — chain A of PDB 4C23 starts at -5."""
    region = _config(viz.build_3dmol_view_html(PDB, [(-5, -2)], chain="A"))["regions"][0]

    assert region["sel"] == {"resi": "-5--2", "chain": "A"}
    assert (region["start"], region["end"]) == (-5, -2)


@pytest.mark.parametrize("options, expected", [
    (dict(crop="all"), None),
    (dict(crop="chain", chain="B"), {"chain": "B"}),
    (dict(crop="chain"), None),                                     # no chain to keep: everything
    (dict(crop="range", crop_range=(30, 10), chain="A"), {"resi": "10-30", "chain": "A"}),
    (dict(crop="epitopes", crop_radius=0),
     {"or": [{"resi": "1-2"}]}),
    (dict(crop="epitopes", crop_radius=6),
     {"within": {"distance": 6.0, "sel": {"or": [{"resi": "1-2"}]}}, "byres": True}),
])
def test_crop_becomes_the_visible_selection(options, expected):
    assert _config(viz.build_3dmol_view_html(PDB, [(1, 2)], **options))["visible"] == expected


def test_cropping_to_the_epitopes_with_none_highlighted_shows_everything():
    """Keeping "only the highlighted regions" of a structure with none would show nothing."""
    assert _config(viz.build_3dmol_view_html(PDB, [], crop="epitopes"))["visible"] is None


def test_labels_start_hidden_when_they_would_cover_the_structure():
    many = [(i, i) for i in range(1, viz.MAX_AUTO_LABELS + 2)]

    assert _config(viz.build_3dmol_view_html(PDB, many[:3]))["labels"] is True
    assert _config(viz.build_3dmol_view_html(PDB, many))["labels"] is False
    # an explicit choice is kept either way
    assert _config(viz.build_3dmol_view_html(PDB, many, show_labels=True))["labels"] is True
    assert _config(viz.build_3dmol_view_html(PDB, many[:3], show_labels=False))["labels"] is False


def test_a_missing_peptide_is_empty_not_the_word_nan():
    df = _df()
    df.loc[0, "Peptide Sequence"] = float("nan")

    assert viz.epitope_list(df, protein="SPIKE", specie="SARS")[0].peptide == ""


def test_a_region_with_no_method_on_record_is_not_colored_as_several_methods():
    regions = viz.color_regions([viz.Region("R1", 1, 5, 1, 5)], "method", {"A": "#123456"})

    assert regions[0].color == "red"


def test_copies_of_the_chain_are_highlighted_too():
    """A dimer colored on one chain only shows half the antigen as if it had no epitopes."""
    alone = _config(viz.build_3dmol_view_html(PDB, [(1, 2)], chain="A"))
    dimer = _config(viz.build_3dmol_view_html(PDB, [(1, 2)], chain="A", copies=["B", "A"]))

    assert alone["regions"][0]["sel"] == {"resi": "1-2", "chain": "A"}
    assert alone["regions"][0]["chains"] == ["A"]
    assert dimer["regions"][0]["sel"] == {"resi": "1-2", "chain": ["A", "B"]}
    assert dimer["regions"][0]["chains"] == ["A", "B"]
    # without a chain every chain is highlighted already, so there is nothing to add
    everywhere = _config(viz.build_3dmol_view_html(PDB, [(1, 2)], copies=["B"]))
    assert everywhere["regions"][0]["chains"] == []
    # cropping to "the mapped chain" still means that one chain
    cropped = _config(viz.build_3dmol_view_html(PDB, [(1, 2)], chain="A", copies=["B"],
                                                crop="chain"))
    assert cropped["visible"] == {"chain": "A"}


def test_a_crop_that_cannot_apply_is_announced_in_the_viewer():
    """Asking for "only the regions" and getting the whole structure, unexplained, looks broken."""
    def note(*args, **kwargs):
        return _config(viz.build_3dmol_view_html(PDB, *args, **kwargs))["note"]

    assert note([(1, 2)]) == ""
    assert note([(1, 2)], crop="epitopes") == ""
    assert "região destacada" in note([], crop="epitopes")
    assert "não se aplica" in note([(1, 2)], crop="chain")          # no chain to keep


def test_the_viewer_is_told_which_format_it_is_reading():
    cif = ("data_X\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n"
           "_atom_site.auth_seq_id\nATOM CA 1\n")

    assert _config(viz.build_3dmol_view_html(PDB, [(1, 2)]))["format"] == "pdb"
    assert _config(viz.build_3dmol_view_html(cif, [(1, 1)]))["format"] == "mcif"
