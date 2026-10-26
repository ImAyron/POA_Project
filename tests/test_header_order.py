"""A inversão `Espécie_Proteína` alcançando as predições, e não só a referência.

O caso real que expôs isto: um consenso de E2 de alphavirus, com cabeçalho
``alphavirus_e2_consenso`` — espécie **alphavirus**, proteína **E2** — e o JSON do BepiPred-2.0
com a mesma string como chave de antígeno.

Marcada a caixa da etapa 1, o ``-f`` era reescrito para proteína ``E2`` / espécie ``ALPHAVIRUS``,
e as predições continuavam a ser lidas ao contrário. As duas metades da mesma execução discordavam,
e a etapa 3 relatava zero coincidência entre elas — depois de mandar marcar a caixa que já estava
marcada.

A correção é num ponto só, depois da consolidação: cada método carrega a convenção à sua maneira
(chave de JSON, cabeçalho FASTA, nome de arquivo), e tratar os seis formatos separadamente seria
seis lugares para manter em sincronia.
"""
import json

import pytest

from poa.core import pipeline
from poa.gui import backend


def _b2_json(path, key="alphavirus_e2_consenso"):
    """Um JSON do BepiPred-2.0 com uma chave de antígeno e um resíduo antigênico."""
    path.write_text(json.dumps({"antigens": {key: {
        "AA": list("MKTAYIA"), "PRED": [0, 0.9, 0.9, 0, 0, 0, 0]}}}), encoding="utf-8")
    return str(path)


def _args(tmp_path, b2, reference_header, swap):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">{reference_header}\nMKTAYIA\n", encoding="utf-8")
    return backend.build_poa1_args({"b2": b2}, {"swap": swap}, str(fasta),
                                   str(tmp_path / "out"))


def test_the_swap_reaches_the_predictions(tmp_path):
    """
    Espécie e proteína saem na ordem certa também nos epítopos.

    Antes, a referência era corrigida e as predições não: o mesmo
    ``alphavirus_e2_consenso`` virava proteína ``E2``/espécie ``ALPHAVIRUS`` de um lado e
    proteína ``ALPHAVIRUS``/espécie ``E2`` do outro.
    """
    b2 = _b2_json(tmp_path / "ep.json")
    args = _args(tmp_path, b2, "E2_ALPHAVIRUS_CONSENSO", swap=True)

    result = pipeline.run_poa1(args)

    # O parser do BepiPred-2.0 preserva a caixa da chave do JSON; o que este teste fixa é a
    # ordem dos campos, não a caixa.
    assert {s.upper() for s in result.predictions["Specie"]} == {"ALPHAVIRUS"}
    assert {s.upper() for s in result.predictions["Protein"]} == {"E2"}


def test_without_the_flag_nothing_is_swapped(tmp_path):
    """O padrão preserva o comportamento da linha de comando, que não tem a caixa."""
    b2 = _b2_json(tmp_path / "ep.json")
    args = _args(tmp_path, b2, "ALPHAVIRUS_E2_CONSENSO", swap=False)

    result = pipeline.run_poa1(args)

    assert {s.upper() for s in result.predictions["Specie"]} == {"E2"}
    assert {s.upper() for s in result.predictions["Protein"]} == {"ALPHAVIRUS"}


def test_the_reference_and_the_predictions_agree_after_the_swap(tmp_path):
    """
    O sintoma que o usuário via: zero coincidência entre o -f e os epítopos.

    É o diagnóstico da etapa 3, e ele só some quando os dois lados são corrigidos.
    """
    b2 = _b2_json(tmp_path / "ep.json")
    args = _args(tmp_path, b2, "E2_ALPHAVIRUS_CONSENSO", swap=True)

    result = pipeline.run_poa1(args)
    diag = backend.species_diagnostics(args.f, result.predictions)

    assert diag["predicted"] == ["ALPHAVIRUS"]
    assert not diag["missing_from_reference"], diag


def test_the_per_species_fasta_uses_the_corrected_species(tmp_path):
    """
    Os arquivos da etapa 4 são nomeados pela espécie, e é por eles que a conservância agrupa.

    Com a ordem trocada o POA1 escrevia ``E2_epitopes.fasta`` — a conservância comparava uma
    "espécie" que era o nome da proteína.
    """
    b2 = _b2_json(tmp_path / "ep.json")
    args = _args(tmp_path, b2, "E2_ALPHAVIRUS_CONSENSO", swap=True)

    result = pipeline.run_poa1(args)
    written = sorted(p.name for p in __import__("pathlib").Path(result.conservancy_dir).glob("*"))

    assert any(name.upper().startswith("ALPHAVIRUS") for name in written), written
    assert not any(name.upper().startswith("E2_") for name in written), written


def test_the_columns_keep_their_order(tmp_path):
    """A troca é de conteúdo, não de esquema: o resto do pipeline lê pelas 7 colunas padrão."""
    from poa.core.models import STANDARD_COLUMNS

    b2 = _b2_json(tmp_path / "ep.json")
    args = _args(tmp_path, b2, "E2_ALPHAVIRUS_CONSENSO", swap=True)

    result = pipeline.run_poa1(args)

    assert list(result.predictions.columns) == list(STANDARD_COLUMNS)


def test_an_empty_result_does_not_break_the_swap(tmp_path):
    """Nenhum epítopo sobrevivendo aos filtros é caso normal, não erro."""
    b2 = _b2_json(tmp_path / "ep.json")
    args = _args(tmp_path, b2, "E2_ALPHAVIRUS_CONSENSO", swap=True)
    args.bmin = 50          # nenhum epítopo chega a 50 resíduos

    result = pipeline.run_poa1(args)

    assert result.predictions.empty


def test_the_cli_exposes_the_flag():
    """A interface e a linha de comando compartilham o core; a opção não pode existir só numa."""
    from poa.cli.poa1 import build_parser

    args = build_parser().parse_args(["-f", "r.fasta", "-d", "out", "-b2", "x.json", "-swap"])
    assert args.swap is True

    default = build_parser().parse_args(["-f", "r.fasta", "-d", "out", "-b2", "x.json"])
    assert default.swap is False


# --------------------------------------------------------------- a prévia da etapa 2
def test_the_step2_preview_follows_the_declared_order(tmp_path):
    """
    A prévia lê o arquivo direto, sem passar pelo POA1 — e precisava ser avisada da inversão.

    Era o sintoma que restava depois de corrigir o pipeline: a etapa 1 concordava com o usuário,
    e a etapa 2 logo adiante ainda anunciava `Espécies: E2` para uma proteína chamada E2. A prévia
    existe para flagrar arquivo errado ou convenção ilegível, e mostrar a ordem que o pipeline
    *não* vai usar trabalha contra isso.
    """
    b2 = _b2_json(tmp_path / "ep.json")

    assert backend.summarize_sources({"b2": [b2]}, swap=False)[0].species == ["E2"]
    assert backend.summarize_sources({"b2": [b2]}, swap=True)[0].species == ["ALPHAVIRUS"]


def test_the_preview_of_a_fasta_follows_the_declared_order(tmp_path):
    """O mesmo vale para os métodos cujo arquivo é FASTA (-b3, -x), não só para o JSON."""
    fasta = tmp_path / "b3.fasta"
    fasta.write_text(">alphavirus_e2_consenso\nmktAYIAkqrq\n", encoding="utf-8")

    plain = backend.summarize_sources({"b3": [str(fasta)]}, swap=False)[0]
    swapped = backend.summarize_sources({"b3": [str(fasta)]}, swap=True)[0]

    assert plain.species == ["E2"] and plain.proteins == ["ALPHAVIRUS"]
    assert swapped.species == ["ALPHAVIRUS"] and swapped.proteins == ["E2"]
    # a tabela por registro acompanha, não só o agregado
    assert swapped.rows[0]["Espécie"] == "ALPHAVIRUS"
    assert swapped.rows[0]["Proteína"] == "E2"
