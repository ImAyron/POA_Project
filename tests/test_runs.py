"""Analysis folders: results/DDMMAAAA-NOME_DO_TESTE, one per test, none overwriting another."""
import json
from datetime import date

from poa.gui import backend


# --------------------------------------------------------------------------- naming
def test_slugify_run_name_handles_accents_and_punctuation():
    assert backend.slugify_run_name("Teste DENV — world set") == "TESTE_DENV_WORLD_SET"
    assert backend.slugify_run_name("conservância 70%") == "CONSERVANCIA_70"
    assert backend.slugify_run_name("  ") == "ANALISE"          # never an empty folder name
    assert len(backend.slugify_run_name("x" * 200)) == 60       # bounded


def test_run_dir_name_puts_the_date_first():
    assert backend.run_dir_name("Teste 1", date(2026, 9, 27)) == "27092026-TESTE_1"


# --------------------------------------------------------------------------- creating
def test_create_run_builds_an_isolated_tree(tmp_path):
    ctx = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))

    assert ctx.root == tmp_path / "27092026-TESTE_A"
    for sub in ("inputs", "mhcii", "poa1_out", "conservancy_csv", "poa2_out"):
        assert (ctx.root / sub).is_dir()

    meta = json.loads((ctx.root / backend.RUN_METADATA_FILENAME).read_text(encoding="utf-8"))
    assert meta["label"] == "Teste A"          # the name as typed, not the slug
    assert meta["created_at"]


def test_two_analyses_do_not_share_files(tmp_path):
    a = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))
    b = backend.create_run(tmp_path, "Teste B", date(2026, 9, 27))
    (a.inputs_dir / "proteins.fasta").write_text(">E_SP_x\nAAA\n")

    assert a.root != b.root
    assert not (b.inputs_dir / "proteins.fasta").exists()


def test_cache_is_shared_across_analyses(tmp_path):
    """Cache keys are content hashes, so a second analysis must not re-submit the same input."""
    a = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))
    b = backend.create_run(tmp_path, "Teste B", date(2026, 9, 27))

    assert a.cache_dir == b.cache_dir == tmp_path / "cache"
    a.cache().put("mhcii", "key1", "cached-tsv", "tsv")
    assert b.cache().get("mhcii", "key1", "tsv") == "cached-tsv"


def test_same_name_same_day_reuses_the_folder(tmp_path):
    first = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))
    (first.inputs_dir / "proteins.fasta").write_text(">E_SP_x\nAAA\n")
    again = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))

    assert again.root == first.root
    assert (again.inputs_dir / "proteins.fasta").exists()   # nothing was wiped


# --------------------------------------------------------------------------- listing
def test_list_runs_is_newest_first(tmp_path):
    backend.create_run(tmp_path, "Antiga", date(2026, 1, 5))
    backend.create_run(tmp_path, "Nova", date(2026, 9, 27))
    backend.create_run(tmp_path, "Meio", date(2026, 5, 10))

    assert [r.label for r in backend.list_runs(tmp_path)] == ["Nova", "Meio", "Antiga"]


def test_list_runs_reports_progress(tmp_path):
    ctx = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))
    (ctx.conservancy_epitopes_dir).mkdir(parents=True, exist_ok=True)
    (ctx.conservancy_epitopes_dir / "SP_epitopes.fasta").write_text(">SP_E_M_1_3\nAYI\n")
    (ctx.conservancy_csv_dir / "SP_conservancy.csv").write_text("Epitope #\n1\n")

    run = backend.list_runs(tmp_path)[0]
    assert (run.n_epitope_fastas, run.n_csvs, run.has_poa2) == (1, 1, False)
    assert "POA1: 1 espécie(s)" in run.progress


def test_legacy_tree_is_listed_last_and_not_moved(tmp_path):
    """The flat results/ written before this convention keeps working exactly where it is."""
    (tmp_path / "poa1_out" / "Conservancy Analysis").mkdir(parents=True)
    (tmp_path / "poa1_out" / "Conservancy Analysis" / "SP_epitopes.fasta").write_text(">a\nAYI\n")
    backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))

    runs = backend.list_runs(tmp_path)
    assert [r.legacy for r in runs] == [False, True]
    legacy = runs[-1]
    assert legacy.path == tmp_path
    assert legacy.n_epitope_fastas == 1
    assert (tmp_path / "poa1_out").is_dir()      # still there, nothing migrated


def test_list_runs_ignores_unrelated_folders(tmp_path):
    backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))
    (tmp_path / "notas").mkdir()
    (tmp_path / "99999999-DATA_INVALIDA").mkdir()

    assert [r.name for r in backend.list_runs(tmp_path)] == ["27092026-TESTE_A"]


def test_list_runs_on_a_missing_root(tmp_path):
    assert backend.list_runs(tmp_path / "nada") == []


# --------------------------------------------------------------------------- reopening
def test_open_run_round_trip(tmp_path):
    created = backend.create_run(tmp_path, "Teste A", date(2026, 9, 27))
    (created.inputs_dir / "proteins.fasta").write_text(">E_SP_x\nAAA\n")

    reopened = backend.open_run(tmp_path, backend.list_runs(tmp_path)[0].path)
    assert reopened.root == created.root
    assert reopened.cache_dir == tmp_path / "cache"
    assert (reopened.inputs_dir / "proteins.fasta").exists()


def test_open_legacy_run_keeps_its_own_cache(tmp_path):
    (tmp_path / "inputs").mkdir()
    ctx = backend.open_run(tmp_path, tmp_path)
    assert ctx.root == tmp_path
    assert ctx.cache_dir == tmp_path / "cache"
