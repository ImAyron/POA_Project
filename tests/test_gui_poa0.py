"""The POA0 search panel inside step 1, driven as a user would, with the databases stubbed.

Step 1 used to have one way in: upload a FASTA you already had. The panel adds the other — name
an organism and let the pipeline find it. Everything here runs offline; the UniProt/PDBe calls
are replaced so the test exercises the form, not the network.
"""
import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from datetime import date  # noqa: E402
from pathlib import Path  # noqa: E402

from streamlit.testing.v1 import AppTest  # noqa: E402

from poa.gui import backend  # noqa: E402
from poa.services import database_client as db  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")

SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 80          # 1600 aa

TAXON = db.Taxon(taxon_id=118655, scientific_name="Oropouche virus")
PROTEIN = db.Protein(accession="A0A0D4BSW3", name="Envelopment polyprotein",
                     length=len(SEQUENCE), reviewed=False, pdb_ids=["6H3X"])
REGION = db.Region(name="Bunyavirus glycoprotein G1", start=514, end=1359, source="Domain+Pfam")
ENTRY = {"sequence": {"length": len(SEQUENCE), "value": SEQUENCE},
         "organism": {"scientificName": "Oropouche virus"},
         "proteinDescription": {"submissionNames": [
             {"fullName": {"value": "Envelopment polyprotein"}}]},
         "features": []}


@pytest.fixture
def stubbed(tmp_path, monkeypatch):
    """An analysis folder plus database answers that never leave the process."""
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    monkeypatch.setattr(db, "search_taxa", lambda name, limit=10: [TAXON])
    monkeypatch.setattr(db, "search_proteins",
                        lambda tid, **kw: [PROTEIN])
    monkeypatch.setattr(db, "get_entry", lambda acc: ENTRY)
    monkeypatch.setattr(db, "list_regions", lambda acc, entry=None, **kw: [REGION])
    monkeypatch.setattr(db, "list_structures", lambda acc, limit=25: [])
    return backend.create_run(tmp_path, "Busca", date(2026, 10, 11))


def _step1(ctx):
    at = AppTest.from_file(APP, default_timeout=90)
    at.session_state["run_path"] = str(ctx.root)
    at.run()
    at.sidebar.radio[0].set_value("1 · Proteínas & Parâmetros").run()
    return at


def _click(at, text):
    next(b for b in at.button if text in b.label).click().run()
    return at


def test_panel_is_offered_without_being_in_the_way(stubbed):
    """
    Someone who already has the FASTA must not have to dismiss a search form first.

    The search is an addition to step 1, not a replacement for it: the panel is collapsed, and
    the uploader, the header-swap checkbox and the POA1 parameters are all still there. The way
    this regressed in practice was not the layout but an exception — anything that raises inside
    the panel aborts the whole step and takes the uploader down with it.
    """
    at = _step1(stubbed)

    assert not at.exception
    panel = next(e for e in at.expander if "buscar pela espécie" in str(e.label).lower())
    assert panel.proto.expanded is False    # closed until asked for
    assert at.file_uploader                 # the original way in, untouched
    assert any("Espécie_Proteína" in str(c.label) for c in at.checkbox)


def test_upload_still_works_when_the_databases_are_unreachable(stubbed, monkeypatch):
    """The old path must not depend on the new one: no network, same step 1."""
    def down(*a, **k):
        raise db.ServiceUnavailable("offline")

    for name in ("search_taxa", "search_proteins", "get_entry", "list_regions", "list_structures"):
        monkeypatch.setattr(db, name, down)

    at = _step1(stubbed)

    assert not at.exception
    assert at.file_uploader


def test_search_narrows_one_question_at_a_time(stubbed):
    at = _step1(stubbed)
    at.text_input(key="poa0_name").set_value("Oropouche virus").run()
    _click(at, "Buscar organismo")

    assert not at.exception
    assert at.session_state["poa0_taxa"] == [TAXON]
    # The protein list is not fetched until the taxon is confirmed: an organism name matches
    # several taxa, so guessing here would query the wrong one.
    assert at.session_state["poa0_proteins"] == []

    _click(at, "Listar proteínas")
    assert at.session_state["poa0_proteins"] == [PROTEIN]
    assert at.session_state["poa0_regions"] == []

    _click(at, "Listar regiões")
    assert at.session_state["poa0_regions"] == [REGION]


def test_generated_fasta_is_adopted_exactly_like_an_upload(stubbed):
    at = _step1(stubbed)
    at.text_input(key="poa0_name").set_value("Oropouche virus").run()
    _click(at, "Buscar organismo")
    _click(at, "Listar proteínas")
    _click(at, "Listar regiões")
    _click(at, "Gerar o FASTA")

    assert not at.exception
    written = Path(at.session_state["ref_raw"][0])
    assert written.is_file()

    lines = written.read_text(encoding="utf-8").splitlines()
    header, sequence = lines[0], "".join(lines[1:])
    # Three fields exactly: Protein_Specie_Accession. More would mean a token smuggled in an
    # underscore, which shifts every field the pipeline reads after it.
    fields = header.lstrip(">").split("_")
    assert len(fields) == 3
    assert fields[1] == "OROPOUCHE"
    assert fields[2] == PROTEIN.accession
    assert len(sequence) in (len(SEQUENCE), REGION.length)

    # The pipeline reads this file the same way it reads an upload.
    rows = backend.parse_fasta_headers(str(written))
    assert rows[0]["ok"] is True
    assert rows[0]["specie"] == "OROPOUCHE"


def test_header_fields_are_editable_before_writing(stubbed):
    """The protein field is what POA2 matches on, so it must not be a name nobody chose."""
    at = _step1(stubbed)
    at.text_input(key="poa0_name").set_value("Oropouche virus").run()
    _click(at, "Buscar organismo")
    _click(at, "Listar proteínas")
    _click(at, "Listar regiões")

    at.text_input(key="poa0_token_name").set_value("GcHead").run()
    at.text_input(key="poa0_token_specie").set_value("OROV").run()
    _click(at, "Gerar o FASTA")

    written = Path(at.session_state["ref_raw"][0])
    assert written.name.startswith("GcHead_OROV")
    assert written.read_text(encoding="utf-8").splitlines()[0].startswith(">GcHead_OROV_")


def test_a_database_outage_is_reported_not_raised(stubbed, monkeypatch):
    """The search is the one step that depends on four external services; it may simply be down."""
    def down(*a, **k):
        raise db.ServiceUnavailable("Could not reach https://rest.uniprot.org")

    monkeypatch.setattr(db, "search_taxa", down)

    at = _step1(stubbed)
    at.text_input(key="poa0_name").set_value("Oropouche virus").run()
    _click(at, "Buscar organismo")

    assert not at.exception                     # the app stays usable
    assert any("falhou" in str(e.value).lower() for e in at.error)
