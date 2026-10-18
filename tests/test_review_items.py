"""Três defeitos silenciosos achados na revisão das etapas 4 e 5.

Nenhum dos três levantava erro: um dava um número que parecia razoável, outro escrevia a
topologia de uma sequência na linha de outra, e o terceiro emitia avisos que não chegavam a
lugar nenhum.
"""
import logging
import warnings
from types import SimpleNamespace

import pytest

from poa.core import topology
from poa.core.parsers import conservancy as cp
from poa.services import conservancy_client as cc


# ======================================================================== 1. resolução do -m
@pytest.mark.parametrize("n, m, esperado", [
    (1, 60, 100.0),      # só 0 ou 100 existem
    (2, 60, 100.0),      # 0, 50, 100 -> 60 vira unanimidade
    (2, 50, 50.0),       # cai exatamente num valor alcançável
    (3, 60, pytest.approx(200 / 3)),   # 0, 33.3, 66.7, 100
    (5, 60, 60.0),       # passos de 20: 60 existe
    (30, 60, 60.0),      # conjunto real do Oropouche
])
def test_effective_match_threshold(n, m, esperado):
    """
    A coluna de conservância é ``casados/n``, então só assume n+1 valores.

    Um ``-m`` entre dois deles se comporta como o de cima — e é assim que "a maioria delas"
    vira "todas elas" sem ninguém perceber.
    """
    assert cp.effective_match_threshold(m, n) == esperado


def test_unknown_set_size_says_nothing():
    assert cp.effective_match_threshold(60, 0) is None


def test_small_set_is_reported_with_the_value_that_takes_effect(tmp_path):
    poa1 = tmp_path / "poa1" / "Conservancy Analysis"
    poa1.mkdir(parents=True)
    (poa1 / "SP_epitopes.fasta").write_text(">SP_P_M_1_3\nAYI\n", encoding="utf-8")
    out = tmp_path / "csvs"
    cc.run_conservancy_for_dir(str(poa1), [("P_SP_a", "GGAYIGG"), ("P_SP_b", "GGAYLGG")],
                               threshold=70, out_dir=str(out), operator=">=")

    avisos = cp.check_match_resolution(str(out), 60)

    assert len(avisos) == 1
    assert "2 proteína(s)" in avisos[0]
    assert "TODAS" in avisos[0]          # -m 60 com n=2 exige unanimidade

    # Com um -m que existe na grade, nada a dizer.
    assert cp.check_match_resolution(str(out), 50) == []


def test_poa2_logs_the_resolution_warning(tmp_path, caplog):
    poa1 = tmp_path / "poa1" / "Conservancy Analysis"
    poa1.mkdir(parents=True)
    (poa1 / "SP_epitopes.fasta").write_text(">SP_P_M_1_3\nAYI\n", encoding="utf-8")
    out = tmp_path / "csvs"
    cc.run_conservancy_for_dir(str(poa1), [("P_SP_a", "GGAYIGG"), ("P_SP_b", "GGAYLGG")],
                               threshold=70, out_dir=str(out), operator=">=")

    with caplog.at_level(logging.WARNING):
        cp.EpitConservAnalysis(70, ">=", 60, 100, 0, str(out))

    assert "TODAS" in caplog.text


# ======================================================================== 2. cabeçalhos repetidos
def test_repeated_protein_headers_are_refused(tmp_path, monkeypatch):
    """
    Dois registros com o mesmo cabeçalho faziam a topologia de um ser escrita na linha do outro.

    O dicionário de anotações guardava o ÚLTIMO; o casamento escolhia o PRIMEIRO. Com
    comprimentos diferentes a fatia da anotação também saía curta, e as três frações
    Portion_* somavam menos de 1 sem erro nenhum.
    """
    import pandas as pd

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">P_SP_ref\nAYIQAYIQ\n>P_SP_ref\nQQQQQQQQQQQQ\n", encoding="utf-8")

    monkeypatch.setattr(topology, "pyTMHMMpredict",
                        lambda path: (["P_SP_REF", "P_SP_REF"], ["oooooooo", "oooooooooooo"]))

    df = pd.DataFrame({"Epitope name": ["SP_P_M_1_4"], "Epitope sequence": ["AYIQ"]})
    args = SimpleNamespace(f=str(fasta))

    with pytest.raises(Exception, match="repeated headers"):
        topology.tmhmmAnalysis(args, df)


def test_distinct_headers_still_work(tmp_path, monkeypatch):
    import pandas as pd

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">P_SP_a\nAYIQAYIQ\n>P_SP_b\nQQQQQQQQ\n", encoding="utf-8")
    monkeypatch.setattr(topology, "pyTMHMMpredict",
                        lambda path: (["P_SP_A", "P_SP_B"], ["oooooooo", "oooooooo"]))

    df = pd.DataFrame({"Epitope name": ["SP_P_M_1_4"], "Epitope sequence": ["AYIQ"]})
    out = topology.tmhmmAnalysis(SimpleNamespace(f=str(fasta)), df)

    assert len(out) == 1
    assert out.iloc[0]["Portion_Outside"] == 1.0


# ======================================================================== 3. avisos que somem
def test_warnings_are_wired_into_the_analysis_log(tmp_path):
    """
    ``topology`` reporta suas ressalvas com ``warnings.warn``, que o Streamlit não mostra e o
    ``poa.log`` não via. Uma ressalva que ninguém consegue ler é o mesmo que nenhuma ressalva.

    O teste confere a ligação, não a entrega: o pytest troca o ``warnings.showwarning`` para
    capturar avisos, o que desfaz o ``captureWarnings`` justamente enquanto os testes rodam.
    O que o código controla é o roteamento — que o ``py.warnings`` escreva nos mesmos handlers
    do ``poa``, inclusive no arquivo da análise.
    """
    from poa.logging_conf import setup_logging, use_log_file

    log = tmp_path / "poa.log"
    setup_logging()
    handler = use_log_file(str(log))

    warnings_logger = logging.getLogger("py.warnings")
    assert handler in warnings_logger.handlers, "o aviso não chegaria ao poa.log da análise"
    assert warnings_logger.propagate is False, "sem isso cada aviso sairia duplicado"

    # E a entrega em si, exercitada pelo caminho que o logging usa, sem depender do
    # showwarning global que o pytest substitui.
    warnings_logger.warning("Epitope SP_P_M_1_4 is not found in the protein(s) FASTA file.")
    handler.flush()
    assert "is not found in the protein" in log.read_text(encoding="utf-8")


def test_a_non_matching_record_does_not_cry_wolf(tmp_path, monkeypatch, recwarn):
    """
    O aviso "not found" disparava uma vez por registro que não continha o epítopo — incluindo o
    caso normal em que um registro seguinte o contém. Avisos assim ensinam a ignorar avisos.

    O que importa (nenhum registro continha) continua sendo avisado, uma vez só.
    """
    import pandas as pd

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">P_SP_a\nQQQQQQQQ\n>P_SP_b\nAYIQAYIQ\n", encoding="utf-8")
    monkeypatch.setattr(topology, "pyTMHMMpredict",
                        lambda path: (["P_SP_A", "P_SP_B"], ["oooooooo", "oooooooo"]))

    df = pd.DataFrame({"Epitope name": ["SP_P_M_1_4"], "Epitope sequence": ["AYIQ"]})
    out = topology.tmhmmAnalysis(SimpleNamespace(f=str(fasta)), df)

    assert out.iloc[0]["Portion_Outside"] == 1.0           # casou com o P_SP_b
    textos = [str(w.message) for w in recwarn]
    assert not [t for t in textos if "not found" in t], textos
