"""Regressões da publicação em lote: lote vazio e retenção do .poa_history."""
from datetime import datetime, timedelta

import pytest

from poa.core import artifacts
from poa.core.artifacts import staged_outputs


def publish(root, text):
    with staged_outputs(root, ["Analysis_report.txt"]) as stage:
        (stage / "Analysis_report.txt").write_text(text, encoding="utf-8")


def generations(root):
    return sorted(p for p in (root / ".poa_history").iterdir() if p.is_dir())


def archived_reports(root):
    return {p.read_text(encoding="utf-8") for p in root.glob(".poa_history/*/Analysis_report.txt")}


def test_a_body_that_writes_nothing_keeps_the_previous_batch(tmp_path):
    """Trava o achado 1: corpo que não gerou arquivo nenhum apagava todo alvo gerenciado.

    O ``all(...)`` do atalho de no-op é vacuamente verdadeiro sobre um lote vazio, a segunda
    metade da condição derruba o atalho e a execução caía no laço que arquiva — removendo os
    resultados da pasta visível e instalando nada no lugar.
    """
    (tmp_path / "Analysis_report.txt").write_text("resultados caros", encoding="utf-8")
    (tmp_path / "reference.fasta").write_text(">ref\nACDE\n", encoding="utf-8")
    (tmp_path / "anotacoes.txt").write_text("arquivo do usuário", encoding="utf-8")

    with staged_outputs(tmp_path, ["Analysis_report.txt", "reference.fasta"]):
        pass  # nenhuma escrita: o corpo desistiu sem produzir o lote

    assert (tmp_path / "Analysis_report.txt").read_text(encoding="utf-8") == "resultados caros"
    assert (tmp_path / "reference.fasta").read_text(encoding="utf-8") == ">ref\nACDE\n"
    assert (tmp_path / "anotacoes.txt").read_text(encoding="utf-8") == "arquivo do usuário"
    assert not (tmp_path / ".poa_history").exists()


def test_an_empty_batch_is_published_only_when_the_caller_says_so(tmp_path):
    """Trava o achado 1: faltava a guarda explícita separando as duas leituras de lote vazio.

    Quem de fato quer publicar um lote vazio — e com isso aposentar o anterior, que é como uma
    reexecução sem espécie alguma limpa os CSVs — passa a declarar isso. Não se deduz mais do
    diretório de stage ter ficado vazio.
    """
    (tmp_path / "Analysis_report.txt").write_text("lote anterior", encoding="utf-8")

    with staged_outputs(tmp_path, ["Analysis_report.txt"], publish_empty=True):
        pass

    assert not (tmp_path / "Analysis_report.txt").exists()
    assert archived_reports(tmp_path) == {"lote anterior"}


def test_history_stops_growing_with_every_rerun(tmp_path):
    """Trava o achado 2: cada publicação acrescentava uma geração e nada podava.

    Dez reexecuções deixavam nove cópias completas do lote sob ``.poa_history``, crescendo
    linearmente com as reexecuções e sem nenhum limite de retenção.
    """
    for i in range(10):
        publish(tmp_path, f"execução {i}")

    kept = generations(tmp_path)
    assert len(kept) <= 3, f"{len(kept)} gerações retidas"
    assert len(kept) == artifacts.HISTORY_GENERATIONS


def test_history_keeps_the_most_recent_generations(tmp_path):
    """Trava o achado 2: a poda tem de descartar as gerações velhas, não as recentes.

    O arquivo só serve para buscar à mão o lote que a reexecução acabou de substituir; geração
    antiga é volume morto, já que nada no projeto restaura a partir dela.
    """
    for i in range(6):
        publish(tmp_path, f"execução {i}")

    assert archived_reports(tmp_path) == {"execução 2", "execução 3", "execução 4"}


def test_a_clock_that_steps_backwards_never_prunes_the_batch_just_replaced(tmp_path, monkeypatch):
    """Trava o achado 2: a poda não pode descartar a única cópia do lote recém-substituído.

    O nome da geração começa pelo instante em que foi criada, então ordenar por nome é ordenar
    por idade — até o relógio de parede andar para trás (NTP, fim de horário de verão) e a
    geração mais nova passar a ordenar como a mais velha.
    """
    class BackwardsClock:
        def __init__(self, start):
            self.value = start

        def now(self):
            self.value -= timedelta(seconds=1)
            return self.value

    monkeypatch.setattr(artifacts, "datetime", BackwardsClock(datetime(2026, 10, 17, 12, 0, 0)))
    for i in range(12):
        publish(tmp_path, f"execução {i}")

    kept = generations(tmp_path)
    assert len(kept) <= 4, f"{len(kept)} gerações retidas"
    assert "execução 10" in archived_reports(tmp_path)


def test_pruning_never_runs_on_the_rollback_path(tmp_path, monkeypatch):
    """Trava o achado 2: a poda não pode rodar quando a publicação falhou.

    O arquivo da geração em curso é de onde o rollback restaura; podar antes de a publicação
    se confirmar trocaria um erro recuperável por perda de dados.
    """
    for i in range(5):
        publish(tmp_path, f"execução {i}")
    before = {p.name for p in generations(tmp_path)}

    replace = artifacts.os.replace
    pruned = []

    def fail_on_install(src, dest):
        if artifacts.Path(src).parent.name.startswith(".poa_stage_"):
            raise OSError("publication failed")
        return replace(src, dest)

    monkeypatch.setattr(artifacts.os, "replace", fail_on_install)
    monkeypatch.setattr(artifacts, "_prune_history", lambda *a, **k: pruned.append(a))
    with pytest.raises(OSError, match="publication failed"):
        publish(tmp_path, "execução 5")

    assert pruned == []
    assert (tmp_path / "Analysis_report.txt").read_text(encoding="utf-8") == "execução 4"
    assert before <= {p.name for p in generations(tmp_path)}
