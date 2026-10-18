"""Defeitos do code-review das correções da auditoria, na camada `poa/gui/app.py`.

Cada caso dirige a página Streamlit como um usuário dirige — estado de sessão, clique de botão,
mensagem na tela — porque os quatro defeitos só aparecem na borda entre o backend, que agora
levanta e recusa mais do que antes, e a página, que precisa continuar de pé e dizer a coisa certa.
"""
import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from pathlib import Path  # noqa: E402
from unittest.mock import Mock  # noqa: E402

from poa.core.parsers import mhcii  # noqa: E402
from poa.gui import backend, viz  # noqa: E402
from poa.services import mhcii_client  # noqa: E402
from poa.services.base import ServiceResult  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")

STEP1 = "1 · Proteínas & Parâmetros"
STEP2 = "2 · Predições (auto/manual)"
STEP7 = "7 · Visualização (2D/3D)"

PEPTIDE = "ACDEFGHIKLMNPQR"


def put(root, name, text):
    path = Path(root) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def tsv(sequence=PEPTIDE):
    return ("allele\tseq_num\tstart\tend\tpeptide\tic50\n"
            f"HLA-DRB1*01:01\t1\t1\t{len(sequence)}\t{sequence}\t10\n")


def pdb(sequence="ACDE", chain="A"):
    names = dict(A="ALA", C="CYS", D="ASP", E="GLU")
    return "\n".join(f"ATOM  {i:>5}  CA  {names[aa]} {chain}{i:>4}    "
                     f"{float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C"
                     for i, aa in enumerate(sequence, 1)) + "\n"


def open_app(ctx, monkeypatch, **state):
    """A página aberta na análise `ctx`, sem log em arquivo e sem asserção de exceção.

    A asserção fica com cada teste: parte deles existe justamente porque a página morria com
    traceback no estado montado aqui.
    """
    from streamlit.testing.v1 import AppTest

    from poa import logging_conf

    monkeypatch.setattr(logging_conf, "use_log_file", lambda _: None)
    at = AppTest.from_file(APP, default_timeout=120)
    at.session_state["run_path"] = str(ctx.root)
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    return at


def record_predictions(monkeypatch):
    """Substitui a chamada ao IEDB e devolve a lista das sequências efetivamente submetidas."""
    calls = []

    def predict(sequence_text, alleles, length=15, method="nn_align", cache=None, **kwargs):
        calls.append(sequence_text)
        return ServiceResult(content=tsv(), source="api")

    monkeypatch.setattr(mhcii_client, "predict", predict)
    return calls


def ten_references(ninth):
    """Dez registros de referência, o nono com o cabeçalho sob teste."""
    headers = [f"E_DENV{i}_ref" for i in range(1, 11)]
    headers[8] = ninth
    return "".join(f">{h}\n{PEPTIDE}\n" for h in headers)


# --------------------------------------------------------------------- achado 1 (app.py:188)
def test_conflicting_reference_fastas_are_reported_instead_of_crashing_the_step(tmp_path, monkeypatch):
    """
    Trava: `merge_fastas` passou a levantar e `_refresh_reference` não pegava nada.

    Dois FASTAs de referência com o mesmo ID e sequências diferentes derrubavam a etapa 1 com
    traceback — e, como os envios ficam no estado da sessão, em TODA reexecução seguinte.
    """
    ctx = backend.WorkContext(tmp_path / "run")
    first = put(ctx.inputs_dir, "ref_0.fasta", ">E_DENV1_ref\nAAAA\n")
    second = put(ctx.inputs_dir, "ref_1.fasta", ">E_DENV1_ref\nCCCC\n")

    at = open_app(ctx, monkeypatch, ref_raw=[first, second])

    assert not at.exception
    assert any("E_DENV1_ref" in str(e.value) for e in at.error)
    assert at.file_uploader                      # a etapa continua utilizável
    assert at.session_state["proteins_path"] is None

    # A segunda reexecução é a que o usuário realmente vê: o estado persiste, a união é refeita.
    at.run()
    assert not at.exception
    assert any("E_DENV1_ref" in str(e.value) for e in at.error)


def test_conflicting_reference_keeps_the_previous_merged_input(tmp_path, monkeypatch):
    """Um envio conflitante não pode apagar o `-f` que já estava valendo."""
    ctx = backend.WorkContext(tmp_path / "run")
    good = put(ctx.inputs_dir, "ref_0.fasta", ">E_DENV1_ref\nAAAA\n")
    clash = put(ctx.inputs_dir, "ref_1.fasta", ">E_DENV1_ref\nCCCC\n")
    merged = put(ctx.inputs_dir, "proteins_all.fasta", ">E_DENV1_ref\nAAAA\n")

    at = open_app(ctx, monkeypatch, ref_raw=[good, clash], proteins_path=merged)

    assert not at.exception
    assert at.session_state["proteins_path"] == merged
    assert Path(merged).read_text(encoding="utf-8") == ">E_DENV1_ref\nAAAA\n"


# --------------------------------------------------------------------- achado 2 (app.py:660)
def test_a_new_mhcii_batch_leaves_behind_nothing_the_reader_would_consolidate(tmp_path, monkeypatch):
    """
    Trava: o conjunto `managed` do staging era `glob("*.html")`, raso e sensível a maiúsculas.

    O leitor (`files_map_MHCIIBD`) percorre subdiretórios e casa `.html` sem distinguir caixa, de
    modo que `sub/OLD2_SP.html` e `OLD_SP.HTML` sobreviviam ao lote novo e entravam no POA1
    junto com ele — a mistura de lotes que o staging existe para impedir.
    """
    ctx = backend.WorkContext(tmp_path / "run")
    refs = put(ctx.inputs_dir, "ref.fasta", f">E_DENV1_ref\n{PEPTIDE}\n")
    put(ctx.mhcii_dir, "OLD_SP.HTML", mhcii_client.normalize_tsv(tsv("AAAAAAAAAAAAAAA")))
    put(ctx.mhcii_dir, "sub/OLD2_SP.html", mhcii_client.normalize_tsv(tsv("CCCCCCCCCCCCCCC")))
    record_predictions(monkeypatch)

    at = open_app(ctx, monkeypatch, proteins_path=refs)
    at.sidebar.radio[0].set_value(STEP2).run()
    next(b for b in at.button if "Automatizar MHC-II" in b.label).click().run()

    assert not at.exception
    species, proteins, paths = mhcii.files_map_MHCIIBD(str(ctx.mhcii_dir))
    assert [Path(p).name for p in paths] == ["E_DENV1_REF.html"]
    assert (species, proteins) == (["DENV1"], ["E"])
    # Arquivados, não apagados: o lote anterior continua recuperável em `.poa_history`.
    archived = {p.name for p in (ctx.mhcii_dir / ".poa_history").rglob("*")
                if p.suffix.lower() == ".html"}
    assert archived == {"OLD_SP.HTML", "OLD2_SP.html"}


# --------------------------------------------------------------------- achado 3 (app.py:1559)
@pytest.mark.parametrize("tamper, expected", [
    ("manifest", "execute o POA1 novamente"),
    ("ambiguous", "referências diferentes"),
])
def test_a_refused_reference_is_not_also_reported_as_missing(tmp_path, monkeypatch, tamper, expected):
    """
    Trava: a etapa 7 mostrava a exceção e, logo abaixo, um aviso que a contradizia.

    `reference_sequence` levanta quando acha a referência e a recusa; o segundo aviso dizia que
    ela não estava em `inputs/` e mandava procurá-la lá — pasta errada, remédio errado.
    """
    ctx = backend.WorkContext(tmp_path / "run")
    refs = put(ctx.inputs_dir, "ref.fasta", ">E_A_ref\nACDE\n")
    pred = put(ctx.inputs_dir, "pred.fasta", ">E_A_M_1_4\nACDE\n")
    result = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    if tamper == "manifest":
        (ctx.poa1_out_dir / "reference.fasta").write_text(">E_A_ref\nAAAA\n", encoding="utf-8")
    else:
        (ctx.poa1_out_dir / "poa1_manifest.json").unlink()
        Path(refs).write_text(">E_A_ref\nACDE\n>E_A_other\nCCCC\n", encoding="utf-8")
    structure = put(ctx.inputs_dir, "structure.pdb", pdb())
    monkeypatch.setattr(viz, "build_3dmol_view_html", Mock(return_value="<div>structure</div>"))

    at = open_app(ctx, monkeypatch, proteins_path=refs, poa1_sources={"x": [pred]},
                  poa1_result=result, pdb_path=structure)
    at.sidebar.radio[0].set_value(STEP7).run()

    assert not at.exception
    messages = [str(m.value) for m in at.warning] + [str(m.value) for m in at.error]
    assert any(expected in m for m in messages)
    assert not any("Não encontrei em" in m for m in messages)


# --------------------------------------------------------------------- achado 4 (app.py:664)
@pytest.mark.parametrize("ninth", ["SEQUENCE", "E_DENV1_ref"])
def test_mhcii_refuses_an_unusable_reference_before_calling_the_iedb(tmp_path, monkeypatch, ninth):
    """
    Trava: os cabeçalhos eram validados DENTRO do laço, depois de a rede já ter sido chamada.

    Um cabeçalho sem espécie, ou repetido, no nono de dez registros custava oito idas ao IEDB
    antes de a entrada ser recusada — e o staging jogava fora as oito. A validação vem antes de
    qualquer chamada, que é o que torna o contrato tudo-ou-nada verdadeiro.
    """
    ctx = backend.WorkContext(tmp_path / "run")
    refs = put(ctx.inputs_dir, "ref.fasta", ten_references(ninth))
    calls = record_predictions(monkeypatch)

    at = open_app(ctx, monkeypatch, proteins_path=refs)
    at.sidebar.radio[0].set_value(STEP2).run()
    next(b for b in at.button if "Automatizar MHC-II" in b.label).click().run()

    assert not at.exception
    assert calls == []
    assert any(ninth in str(e.value) for e in at.error)
    assert list(ctx.mhcii_dir.rglob("*.html")) == []


def test_mhcii_still_accepts_the_two_field_header_the_rest_of_the_pipeline_accepts(tmp_path, monkeypatch):
    """
    Trava: a validação nova exigia três campos onde a etapa 1 e o leitor aceitam dois.

    `Proteína_Espécie` é o nome que a própria etapa 2 pede no envio manual e que o POA1 lê com
    `ID_Sequence` igual a `-`; recusá-lo aqui bloqueava a automação de uma entrada legítima.
    """
    ctx = backend.WorkContext(tmp_path / "run")
    refs = put(ctx.inputs_dir, "ref.fasta", ten_references("E_DENV1"))
    calls = record_predictions(monkeypatch)

    at = open_app(ctx, monkeypatch, proteins_path=refs)
    at.sidebar.radio[0].set_value(STEP2).run()
    next(b for b in at.button if "Automatizar MHC-II" in b.label).click().run()

    assert not at.exception
    assert len(calls) == 10
    assert (ctx.mhcii_dir / "E_DENV1.html").is_file()
    species, proteins, paths = mhcii.files_map_MHCIIBD(str(ctx.mhcii_dir))
    assert len(paths) == 10
    assert species.count("DENV1") == 2        # o de dois campos e o `E_DENV1_ref`
