"""Apagar uma análise: o que some, e o que o código se recusa a apagar.

A exclusão é irreversível e não passa pela lixeira, então as salvaguardas valem mais do que a
funcionalidade. A confirmação na interface (digitar o nome da pasta) é a primeira barreira; estes
testes cobrem a segunda, que é o próprio backend recusar alvos errados mesmo quando chamado com
os argumentos errados.
"""
from datetime import date

import pytest

from poa.gui import backend


@pytest.fixture
def root(tmp_path):
    return tmp_path


def _run_with_files(root, nome, quando=date(2026, 10, 12)):
    ctx = backend.create_run(root, nome, quando)
    (ctx.inputs_dir / "proteins.fasta").write_text(">P_SP_ref\nAYIQ\n", encoding="utf-8")
    (ctx.conservancy_csv_dir / "SP_conservancy.csv").write_text("Epitope #\n1\n", encoding="utf-8")
    (ctx.root / "poa.log").write_text("log\n", encoding="utf-8")
    return ctx


# --------------------------------------------------------------------------- o caminho normal
def test_deleting_removes_the_folder_and_reports_what_went(root):
    ctx = _run_with_files(root, "Teste A")
    assert ctx.root.is_dir()

    removed = backend.delete_run(root, ctx.root)

    assert not ctx.root.exists()
    assert removed["folder"] == ctx.root.name
    assert removed["files"] >= 4                 # run.json + os três que escrevemos
    assert removed["bytes"] > 0


def test_other_analyses_are_untouched(root):
    keep = _run_with_files(root, "Fica")
    drop = _run_with_files(root, "Sai")

    backend.delete_run(root, drop.root)

    assert keep.root.is_dir()
    assert (keep.inputs_dir / "proteins.fasta").is_file()
    assert [r.label for r in backend.list_runs(root)] == ["Fica"]


# --------------------------------------------------------------------------- as recusas
def test_the_results_root_itself_is_refused(root):
    _run_with_files(root, "Teste A")

    with pytest.raises(ValueError, match="raiz"):
        backend.delete_run(root, root)
    assert root.is_dir()


def test_a_path_outside_the_root_is_refused(root, tmp_path):
    """Um argumento errado não pode sair da árvore de resultados."""
    fora = tmp_path.parent / "fora_do_root"
    fora.mkdir(exist_ok=True)
    (fora / "importante.txt").write_text("nao apague", encoding="utf-8")

    with pytest.raises(ValueError, match="fora de"):
        backend.delete_run(root, fora)
    assert (fora / "importante.txt").is_file()


def test_the_shared_cache_is_refused(root):
    """O cache é de todas as análises; apagá-lo junto com uma delas custaria horas de API."""
    ctx = _run_with_files(root, "Teste A")
    cache = root / "cache"
    cache.mkdir(exist_ok=True)
    (cache / "mhcii.tsv").write_text("cached", encoding="utf-8")

    with pytest.raises(ValueError, match="cache"):
        backend.delete_run(root, cache)
    assert (cache / "mhcii.tsv").is_file()
    assert ctx.root.is_dir()


def test_a_folder_that_is_not_an_analysis_is_refused(root):
    """
    ``results/`` também guarda saídas do POA0 e pastas de rascunho; nenhuma delas é uma análise.

    Sem este teste, selecionar a pasta errada apagaria dados que o pipeline nunca criou.
    """
    qualquer = root / "teste_poa0"
    qualquer.mkdir()
    (qualquer / "GcHead_OROV.fasta").write_text(">GcHead_OROV_X\nAYIQ\n", encoding="utf-8")

    with pytest.raises(ValueError, match="não parece uma análise"):
        backend.delete_run(root, qualquer)
    assert (qualquer / "GcHead_OROV.fasta").is_file()


def test_a_missing_folder_says_so(root):
    with pytest.raises(FileNotFoundError):
        backend.delete_run(root, root / "12102026-NAO_EXISTE")


def test_a_file_is_not_a_folder(root):
    arquivo = root / "arquivo.txt"
    arquivo.write_text("x", encoding="utf-8")

    with pytest.raises(ValueError, match="não é uma pasta"):
        backend.delete_run(root, arquivo)
    assert arquivo.is_file()


# --------------------------------------------------------------------------- formatação
@pytest.mark.parametrize("nbytes, esperado", [
    (0, "0 B"), (512, "512 B"), (2048, "2.0 KB"), (5 * 1024 ** 2, "5.0 MB"),
])
def test_human_size(nbytes, esperado):
    assert backend.human_size(nbytes) == esperado


# --------------------------------------------------------------------------- pela interface
pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from pathlib import Path  # noqa: E402

from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")


def _picker(root, monkeypatch):
    """A tela de escolha de análise, com o botão 'Abrir existente' selecionado."""
    monkeypatch.setenv("POA_RESULTS_DIR", str(root))
    at = AppTest.from_file(APP, default_timeout=90)
    at.run()
    return at


def test_delete_is_offered_but_closed(root, monkeypatch):
    _run_with_files(root, "Teste A")
    at = _picker(root, monkeypatch)

    painel = next(e for e in at.expander if "apagar" in str(e.label).lower())
    assert painel.proto.expanded is False       # não é a primeira coisa que se vê


def test_the_button_is_dead_until_the_folder_name_is_typed(root, monkeypatch):
    """
    A confirmação é digitar o nome, não um segundo clique.

    Um par de botões vira reflexo; digitar o nome é uma decisão, e torna impossível apagar uma
    análise diferente daquela que se está lendo na tela.
    """
    ctx = _run_with_files(root, "Teste A")
    at = _picker(root, monkeypatch)

    botao = next(b for b in at.button if "Apagar definitivamente" in b.label)
    assert botao.disabled is True

    at.text_input(key=f"del_{ctx.root.name}").set_value("qualquer coisa").run()
    botao = next(b for b in at.button if "Apagar definitivamente" in b.label)
    assert botao.disabled is True               # nome errado não serve
    assert ctx.root.is_dir()

    at.text_input(key=f"del_{ctx.root.name}").set_value(ctx.root.name).run()
    botao = next(b for b in at.button if "Apagar definitivamente" in b.label)
    assert botao.disabled is False


def test_confirming_deletes_and_returns_to_the_picker(root, monkeypatch):
    a = _run_with_files(root, "Teste A")
    b = _run_with_files(root, "Teste B")
    at = _picker(root, monkeypatch)

    # Qual das duas o seletor abre depende da ordenação (as duas são do mesmo dia), e o painel
    # pertence à análise selecionada — é esse acoplamento que impede apagar a errada.
    campo = next(e for e in at.text_input if str(e.key or "").startswith("del_"))
    alvo = Path(campo.key[len("del_"):])
    alvo = a.root if alvo.name == a.root.name else b.root
    sobrevivente = b.root if alvo == a.root else a.root

    campo.set_value(alvo.name).run()
    next(bt for bt in at.button if "Apagar definitivamente" in bt.label).click().run()

    assert not at.exception
    assert not alvo.exists()
    assert sobrevivente.is_dir()                # a outra análise continua lá
    assert at.session_state.get("run_path") is None
