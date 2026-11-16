"""Step 7 holds several structures: listing them, fetching one by PDB id, choosing among them.

The step used to know a single structure — the last one to arrive — so looking at the same
epitopes on a second PDB entry meant uploading over the first, and uploading the first again to
go back. The files were always kept; what was missing was treating them as a list to choose from.
"""
from pathlib import Path
from unittest.mock import Mock

import pytest

from poa.gui import backend, viz
from poa.services import database_client
from poa.services.base import ServiceUnavailable

STEP7 = "7 · Visualização (2D/3D)"
THREE = dict(A="ALA", C="CYS", D="ASP", E="GLU")


def pdb(sequence="ACDE", chain="A", first=1):
    return "\n".join(
        f"ATOM  {i:>5}  CA  {THREE[aa]} {chain}{first + i - 1:>4}    "
        f"{float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C"
        for i, aa in enumerate(sequence, 1)) + "\n"


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    """
    No test here reaches the network, whatever it forgets to replace.

    ``backend.fetch_structure`` falls back to mmCIF when the PDB-format download fails, so a test
    that fails the first download on purpose would otherwise go and fetch the second for real.
    """
    monkeypatch.setattr(database_client, "_get_text",
                        Mock(side_effect=ServiceUnavailable("sem rede nos testes")))
    monkeypatch.setattr(database_client, "_get_json",
                        Mock(side_effect=ServiceUnavailable("sem rede nos testes")))


@pytest.fixture
def ctx(tmp_path):
    return backend.WorkContext(tmp_path / "run")


def _save(ctx, name, text, **provenance):
    return backend.save_upload(ctx.inputs_dir, name, text.encode("utf-8"), ctx=ctx, **provenance)


# --------------------------------------------------------------------------- the list
def test_every_structure_of_the_analysis_is_listed_oldest_first(ctx):
    first = _save(ctx, "structure_6H3X.pdb", pdb(), role="structure")
    _save(ctx, "ref_0_x.fasta", ">E_A_ref\nACDE\n", role="reference")
    second = _save(ctx, "structure_5U0B.pdb", pdb(), role="structure")

    assert backend.structure_files(ctx) == [first, second]


def test_a_structure_that_predates_the_record_is_listed_too(ctx):
    """Files already on disk have no history; the name POA0 and the upload give them is enough."""
    recorded = _save(ctx, "structure_NOVA.pdb", pdb(), role="structure")
    legacy = ctx.inputs_dir / "structure_ANTIGA.pdb"
    legacy.write_text(pdb(), encoding="utf-8")
    (ctx.inputs_dir / "notas.txt").write_text("nada", encoding="utf-8")

    assert backend.structure_files(ctx) == [recorded, str(legacy)]


def test_a_deleted_structure_leaves_the_list(ctx):
    kept = _save(ctx, "structure_A.pdb", pdb(), role="structure")
    gone = _save(ctx, "structure_B.pdb", pdb(), role="structure")
    Path(gone).unlink()

    assert backend.structure_files(ctx) == [kept]
    assert backend.structure_files(backend.WorkContext(ctx.root.parent / "vazia")) == []


# --------------------------------------------------------------------------- fetch by id
def test_fetching_by_id_saves_and_records_the_structure(ctx, monkeypatch):
    asked = []
    monkeypatch.setattr(database_client, "fetch_structure",
                        lambda code: asked.append(code) or pdb())

    path = backend.fetch_structure(ctx, " 6h3x ")

    assert asked == ["6H3X"]
    assert Path(path).name == "structure_6H3X.pdb"            # the name POA0 uses for the same file
    assert Path(path).read_text(encoding="utf-8") == pdb()
    assert backend.structure_files(ctx) == [path]
    assert backend.restore_analysis(ctx).structure == path    # and it comes back on reopening


def test_an_entry_already_in_the_analysis_is_not_downloaded_again(ctx, monkeypatch):
    """Step 7 must keep working where the first download was the last moment with a network."""
    monkeypatch.setattr(database_client, "fetch_structure", lambda code: pdb())
    first = backend.fetch_structure(ctx, "6H3X")
    monkeypatch.setattr(database_client, "fetch_structure",
                        Mock(side_effect=ServiceUnavailable("offline")))

    assert backend.fetch_structure(ctx, "6h3x") == first


@pytest.mark.parametrize("bad", ["", "6H3", "ABCD", "6H3XX", "../x", "6H 3"])
def test_a_malformed_id_is_refused_before_any_request(ctx, monkeypatch, bad):
    fetch = Mock()
    monkeypatch.setattr(database_client, "fetch_structure", fetch)

    with pytest.raises(ValueError, match="identificador do PDB"):
        backend.fetch_structure(ctx, bad)
    fetch.assert_not_called()
    assert backend.structure_files(ctx) == []


def test_a_failed_download_leaves_nothing_behind(ctx, monkeypatch):
    monkeypatch.setattr(database_client, "fetch_structure",
                        Mock(side_effect=ServiceUnavailable("404")))

    with pytest.raises(ServiceUnavailable):
        backend.fetch_structure(ctx, "9ZZZ")
    assert backend.structure_files(ctx) == []
    assert backend.read_input_history(ctx) == []


# --------------------------------------------------------------------------- the choice
def test_the_chosen_structure_is_what_reopening_restores(ctx):
    """
    Arrivals say which file came last; only the step state says which one was being looked at.

    Going back to an earlier structure records no arrival, so without the choice the analysis
    reopened on the structure nobody had on screen.
    """
    earlier = _save(ctx, "structure_6H3X.pdb", pdb(), role="structure")
    later = _save(ctx, "structure_5U0B.pdb", pdb(), role="structure")
    assert backend.restore_analysis(ctx).structure == later

    backend.record_step_state(ctx, structure=earlier)

    assert backend.restore_analysis(ctx).structure == earlier


def test_a_chosen_structure_that_is_gone_falls_back_to_the_arrivals(ctx):
    earlier = _save(ctx, "structure_6H3X.pdb", pdb(), role="structure")
    later = _save(ctx, "structure_5U0B.pdb", pdb(), role="structure")
    backend.record_step_state(ctx, structure=earlier)
    Path(earlier).unlink()

    assert backend.restore_analysis(ctx).structure == later


# --------------------------------------------------------------------------- the interface
def _open_step7(ctx, monkeypatch, **extra):
    pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    from poa import logging_conf

    monkeypatch.setattr(logging_conf, "use_log_file", lambda _: None)
    refs = ctx.inputs_dir / "ref.fasta"
    refs.parent.mkdir(parents=True, exist_ok=True)
    refs.write_text(">E_A_ref\nACDE\n", encoding="utf-8")
    pred = ctx.inputs_dir / "pred.fasta"
    pred.write_text(">E_A_M_1_4\nACDE\n", encoding="utf-8")
    result = backend.run_poa1(ctx, {"x": str(pred)}, {}, str(refs))

    render = Mock(return_value="<div>structure</div>")
    monkeypatch.setattr(viz, "build_3dmol_view_html", render)
    at = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "poa/gui/app.py"),
                           default_timeout=90)
    state = {"run_path": str(ctx.root), "proteins_path": str(refs),
             "poa1_sources": {"x": [str(pred)]}, "poa1_result": result}
    state.update(extra)
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    at.sidebar.radio[0].set_value(STEP7).run()
    assert not at.exception
    return at, render


def _select(at, label):
    return next(s for s in at.selectbox if s.label == label)


def test_step7_offers_every_structure_and_switches_between_them(ctx, monkeypatch):
    """Two entries of the same protein, numbered differently: the same epitope on each."""
    first = _save(ctx, "structure_AAAA.pdb", pdb(first=1), role="structure")
    second = _save(ctx, "structure_BBBB.pdb", pdb(first=482), role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=second)

    picker = _select(at, "Estrutura em uso")
    assert list(picker.options) == ["AAAA.pdb", "BBBB.pdb"]
    assert picker.value == second
    assert render.call_args.args[1] == [(482, 485)]           # remapped to that file's numbering

    at = picker.set_value(first).run()

    assert not at.exception
    assert at.session_state["pdb_path"] == first
    assert render.call_args.args[1] == [(1, 4)]
    assert backend.restore_analysis(ctx).structure == first   # the choice outlives the session


def test_step7_names_the_regions_it_highlights(ctx, monkeypatch):
    path = _save(ctx, "structure_AAAA.pdb", pdb(first=482), role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)

    regions = render.call_args.kwargs["regions"]
    assert [(r.name, r.seq_start, r.seq_end, r.pdb_start, r.pdb_end) for r in regions] == \
        [("R1", 1, 4, 482, 485)]
    assert render.call_args.kwargs["chain"] == "A"
    table = next(d.value for d in at.dataframe if "Região" in d.value.columns)
    assert list(table["Região"]) == ["R1"] and list(table["PDB"]) == ["A:482–485"]


def test_step7_passes_the_crop_and_the_display_options_to_the_viewer(ctx, monkeypatch):
    path = _save(ctx, "structure_AAAA.pdb", pdb() + pdb("CCCC", "B"), role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)
    assert render.call_args.kwargs["crop"] == "all"

    at = _select(at, "Recorte — o que mostrar da estrutura").set_value("range").run()
    at = next(n for n in at.number_input if n.label == "Do resíduo").set_value(2).run()
    at = next(n for n in at.number_input if n.label == "Até o resíduo").set_value(3).run()
    at = _select(at, "Colorir por").set_value("single").run()
    at = next(c for c in at.checkbox if c.label == "Rótulos nas regiões").set_value(False).run()

    assert not at.exception
    options = render.call_args.kwargs
    assert options["crop"] == "range" and options["crop_range"] == (2, 3)
    assert options["show_labels"] is False
    assert [r.color for r in options["regions"]] == ["red"]
    # and the form survives leaving the step, like every other one
    at = at.sidebar.radio[0].set_value("6 · Resultados").run()
    at = at.sidebar.radio[0].set_value(STEP7).run()
    assert _select(at, "Recorte — o que mostrar da estrutura").value == "range"
    assert render.call_args.kwargs["crop_range"] == (2, 3)


def test_step7_downloads_a_structure_by_its_pdb_id(ctx, monkeypatch):
    monkeypatch.setattr(database_client, "fetch_structure", lambda code: pdb(first=482))
    at, render = _open_step7(ctx, monkeypatch)
    assert render.call_count == 0                             # nothing to show yet

    at = next(t for t in at.text_input if "Identificador do PDB" in t.label).set_value("6h3x").run()
    at = next(b for b in at.button if "Baixar do RCSB" in b.label).click().run()

    assert not at.exception
    assert Path(at.session_state["pdb_path"]).name == "structure_6H3X.pdb"
    assert render.call_args.args[1] == [(482, 485)]


def test_step7_reports_a_failed_download_and_keeps_the_current_structure(ctx, monkeypatch):
    path = _save(ctx, "structure_AAAA.pdb", pdb(), role="structure")
    monkeypatch.setattr(database_client, "fetch_structure",
                        Mock(side_effect=ServiceUnavailable("404 Not Found")))
    at, _ = _open_step7(ctx, monkeypatch, pdb_path=path)

    at = next(t for t in at.text_input if "Identificador do PDB" in t.label).set_value("9ZZZ").run()
    at = next(b for b in at.button if "Baixar do RCSB" in b.label).click().run()

    assert not at.exception
    assert any("9ZZZ" in str(e.value) for e in at.error)
    assert at.session_state["pdb_path"] == path


def test_a_reply_that_is_not_a_structure_is_not_saved_as_one(ctx, monkeypatch):
    """
    A proxy or captive portal answers any URL with a page of its own, status 200.

    An entry already in the analysis is reused instead of downloaded, so a saved error page would
    be "the structure" for every later request — an empty viewer with nothing to explain it.
    """
    monkeypatch.setattr(database_client, "fetch_structure", lambda code: "<html>login</html>")

    with pytest.raises(ServiceUnavailable, match="não é um arquivo PDB"):
        backend.fetch_structure(ctx, "6H3X")
    assert list(ctx.inputs_dir.glob("*")) == []
    assert backend.structure_files(ctx) == []


def test_step7_finds_the_antigen_in_a_complex_and_colors_its_copies(ctx, monkeypatch):
    """Antibody chain first, then the antigen twice: chosen without being asked, both colored."""
    text = pdb("CCCC", "H") + pdb("ACDE", "A", first=482) + pdb("ACDE", "B", first=482)
    path = _save(ctx, "structure_CPLX.pdb", text, role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)

    assert render.call_args.kwargs["chain"] == "A"
    assert render.call_args.kwargs["copies"] == ["B"]
    assert render.call_args.args[1] == [(482, 485)]
    assert any("cadeia(s) B" in str(c.value) for c in at.caption)

    # naming a chain means that chain and nothing else
    at = _select(at, "Cadeia").set_value("B").run()
    assert render.call_args.kwargs["chain"] == "B"
    assert render.call_args.kwargs["copies"] == []

    # and naming the wrong one is still refused, not silently corrected
    at = _select(at, "Cadeia").set_value("H").run()
    assert not at.exception
    assert any("Estrutura incompatível" in str(e.value) for e in at.error)
    assert render.call_args.args[1] == []


def test_step7_crop_range_is_stated_in_the_numbering_of_the_mapped_chain(ctx, monkeypatch):
    """The first chain is numbered 1-4 and is not the protein; the range must offer 482-485."""
    text = pdb("CCCC", "H") + pdb("ACDE", "A", first=482)
    path = _save(ctx, "structure_CPLX.pdb", text, role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)

    at = _select(at, "Recorte — o que mostrar da estrutura").set_value("range").run()

    assert not at.exception
    first = next(n for n in at.number_input if n.label == "Do resíduo")
    assert (first.min, first.max, first.value) == (482, 485, 482)
    assert render.call_args.kwargs["crop_range"] == (482, 485)


def test_step7_does_not_call_an_empty_highlight_a_success(ctx, monkeypatch):
    """Nothing selected to highlight: the plain structure is not announced in green."""
    path = _save(ctx, "structure_PART.pdb", pdb("ACDE", "A"), role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)
    assert any("destacad" in str(s.value) for s in at.success)

    at = next(m for m in at.multiselect if "Métodos" in m.label).set_value([]).run()

    assert not at.exception
    assert render.call_args.args[1] == []
    assert not any("destacad" in str(s.value) for s in at.success)
    assert any("Nenhum método selecionado" in str(i.value) for i in at.info)


def test_step7_fields_keep_every_change_not_every_other_one(ctx, monkeypatch):
    """
    Changing the same field twice in a row used to drop the second change.

    Found on this form in a browser: tick "Mostrar superfície" and it ticks, untick it and it
    stays ticked. See ``_field`` in ``app.py`` for the mechanism; this drives each kind of widget
    the step has through several consecutive changes.
    """
    text = pdb("ACDE", "A") + pdb("ACDE", "B") + pdb("CCCC", "H")
    first = _save(ctx, "structure_AAAA.pdb", text, role="structure")
    second = _save(ctx, "structure_BBBB.pdb", text, role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=first)

    for wanted in ("method", "single", "region", "single"):
        at = _select(at, "Colorir por").set_value(wanted).run()
        assert _select(at, "Colorir por").value == wanted
        assert at.session_state["params_viz"]["color_by"] == wanted
    for wanted in ("B", "H", "A", ""):
        at = _select(at, "Cadeia").set_value(wanted).run()
        assert at.session_state["params_viz"]["chain"] == wanted
    for wanted in (True, False, True, False):
        at = next(c for c in at.checkbox if c.label == "Mostrar superfície").set_value(wanted).run()
        assert render.call_args.kwargs["show_surface"] is wanted
    for wanted in (second, first, second):
        at = _select(at, "Estrutura em uso").set_value(wanted).run()
        assert at.session_state["pdb_path"] == wanted
    assert not at.exception


# --------------------------------------------------------------------------- mmCIF
CIF = ("data_TEST\n#\n_struct.title 'A test structure'\n_exptl.method 'X-RAY DIFFRACTION'\n"
       "_refine.ls_d_res_high 2.10\n#\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n"
       "_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_seq_id\n"
       "_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n"
       "ATOM CA . ALA 1 482 A 1\nATOM CA . CYS 2 483 A 1\nATOM CA . ASP 3 484 A 1\n"
       "ATOM CA . GLU 4 485 A 1\n#\n")


def test_an_entry_with_no_pdb_format_file_is_fetched_as_mmcif(ctx, monkeypatch):
    """
    The largest entries exist only as mmCIF.

    The failed download used to end with "baixe do site e envie o arquivo" — advice that could
    not be followed, because the uploader took ``.pdb`` alone.
    """
    monkeypatch.setattr(database_client, "fetch_structure",
                        Mock(side_effect=ServiceUnavailable("404")))
    monkeypatch.setattr(database_client, "fetch_structure_cif", lambda code: CIF)

    path = backend.fetch_structure(ctx, "9ZZZ")

    assert Path(path).name == "structure_9ZZZ.cif"
    assert backend.structure_files(ctx) == [path]
    # and it is found again without a second download, whichever format it was saved in
    monkeypatch.setattr(database_client, "fetch_structure_cif", Mock(side_effect=AssertionError))
    assert backend.fetch_structure(ctx, "9zzz") == path


def test_the_header_of_an_mmcif_file_says_what_it_is():
    summary = database_client.structure_summary(CIF)

    assert summary["title"] == "A TEST STRUCTURE"
    assert summary["method"] == "X-RAY DIFFRACTION"
    assert summary["resolution"].startswith("2.10")
    assert summary["chains"] == "A"


def test_the_header_is_not_cut_at_a_mention_of_the_coordinate_table():
    """A real entry's revision history lists the ``_atom_site`` items it changed, near the top."""
    mentioned = CIF.replace("#\n_struct.title",
                            "#\n_pdbx_audit_revision_item.item '_atom_site.auth_asym_id'\n"
                            "_struct.title", 1)

    assert database_client.structure_summary(mentioned)["title"] == "A TEST STRUCTURE"


def test_step7_shows_an_mmcif_structure_like_any_other(ctx, monkeypatch):
    path = _save(ctx, "structure_TEST.cif", CIF, role="structure")
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)

    assert render.call_args.args[1] == [(482, 485)]
    assert render.call_args.kwargs["chain"] == "A"
    assert any("A TEST STRUCTURE" in str(i.value) for i in at.info)


# --------------------------------------------------------------------------- SIFTS
def _poa0_report(ctx, structure, *, offset, length=4, chain="A"):
    import json

    report = {"arquivos": {"estrutura": structure}, "sequencia": {"tamanho": length},
              "regiao": {"start": 482, "end": 485},
              "estrutura": {"pdb_id": "AAAA", "offset_regiao": offset, "cobre_regiao": True,
                            "sifts": {"chain_id": chain}}}
    (ctx.inputs_dir / backend.POA0_REPORT_FILENAME).write_text(json.dumps(report), encoding="utf-8")


def test_the_sifts_record_names_its_chain(ctx):
    """The report stores ``chain_id``; ``chain`` was read, and was always empty."""
    path = _save(ctx, "structure_AAAA.pdb", pdb(first=482), role="structure")
    _poa0_report(ctx, path, offset=481, chain="B")

    assert backend.structure_mapping(ctx, path)["chain"] == "B"


def test_step7_uses_the_sifts_mapping_poa0_recorded(ctx, monkeypatch):
    """The antigen is the second chain; SIFTS says which one, and the copies still follow."""
    text = pdb("CCCC", "H") + pdb("ACDE", "A", first=482) + pdb("ACDE", "B", first=482)
    path = _save(ctx, "structure_AAAA.pdb", text, role="structure")
    _poa0_report(ctx, path, offset=481)
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)

    assert render.call_args.args[1] == [(482, 485)]
    assert render.call_args.kwargs["chain"] == "A"
    assert render.call_args.kwargs["copies"] == ["B"]
    assert any("SIFTS" in str(c.value) for c in at.caption)
    assert any("via mapeamento SIFTS" in str(s.value) for s in at.success)


def test_step7_falls_back_to_the_alignment_when_the_sifts_record_does_not_apply(ctx, monkeypatch):
    """
    The offset POA0 recorded is for the region POA0 wrote.

    Against another reference it places nothing, and that says the record does not apply — not
    that the structure is the wrong protein. The alignment answers before anything is refused.
    """
    path = _save(ctx, "structure_AAAA.pdb", pdb(first=482), role="structure")
    _poa0_report(ctx, path, offset=100)                       # some other region's offset
    at, render = _open_step7(ctx, monkeypatch, pdb_path=path)

    assert not at.error
    assert render.call_args.args[1] == [(482, 485)]
    assert not any("SIFTS" in str(c.value) for c in at.caption)
    assert any("via alinhamento" in str(s.value) for s in at.success)


# --------------------------------------------------------------------------- the PDBe search
def _with_accession(ctx):
    refs = ctx.inputs_dir / "ref.fasta"
    refs.parent.mkdir(parents=True, exist_ok=True)
    refs.write_text(">E_A_P29990\nACDE\n", encoding="utf-8")


def test_step7_search_lists_and_downloads_into_the_structure_list(ctx, monkeypatch):
    """The search is a third way in; what it downloads is a structure like the others."""
    found = [database_client.Structure(pdb_id="6H3X", accession="P29990", chains=["A"],
                                       unp_start=1, unp_end=4, coverage=1.0, resolution=2.1)]
    asked = []
    monkeypatch.setattr(database_client, "list_structures",
                        lambda accession, limit=15: asked.append(accession) or found)
    monkeypatch.setattr(database_client, "fetch_structure", lambda code: pdb(first=482))
    at, render = _open_step7(ctx, monkeypatch, proteins_path=None)
    # the reference _open_step7 wrote has no accession; give this analysis one that does
    at.session_state["proteins_path"] = str(ctx.inputs_dir / "ref_acc.fasta")
    (ctx.inputs_dir / "ref_acc.fasta").write_text(">E_A_P29990\nACDE\n", encoding="utf-8")
    at = at.run()

    at = next(b for b in at.button if b.label == "Buscar estruturas").click().run()
    assert asked == ["P29990"]
    at = next(b for b in at.button if "Baixar e usar" in b.label).click().run()

    assert not at.exception
    assert Path(at.session_state["pdb_path"]).name == "structure_6H3X.pdb"
    assert backend.structure_files(ctx) == [at.session_state["pdb_path"]]
    assert render.call_args.args[1] == [(482, 485)]


def test_a_listing_is_not_offered_under_another_proteins_name(ctx, monkeypatch):
    """Search results answer one accession; kept after the protein changed, they mislabel it."""
    at, _ = _open_step7(ctx, monkeypatch)
    at.session_state["proteins_path"] = str(ctx.inputs_dir / "ref_acc.fasta")
    (ctx.inputs_dir / "ref_acc.fasta").write_text(">E_A_P29990\nACDE\n", encoding="utf-8")
    at.session_state["viz_structures"] = [database_client.Structure(pdb_id="4C23", accession="Q00000")]
    at.session_state["viz_structures_for"] = "Q00000"
    at = at.run()

    assert not at.exception
    assert not [b for b in at.button if "Baixar e usar" in b.label]
