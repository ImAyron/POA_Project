"""Defeitos encontrados no code-review das correções da auditoria — área da conservância.

Cada caso trava um defeito que as próprias correções introduziram em
``poa/core/parsers/conservancy.py`` e ``poa/services/conservancy_client.py``, não um defeito da
base original: sem a correção do review, todos falham.
"""
import json

import pytest

from poa.core.parsers import conservancy as cp
from poa.services import conservancy_client as cc

#: Conjunto de comparação com duas proteínas, para o percentual não ser 0% ou 100% por construção.
PROTEINS = [("E_A_ref", "ACDEFGHI"), ("E_A_two", "ACDEFGHL")]


def conservancy_csv(path, threshold=70, operator=">="):
    """Um resultado de conservância de verdade, com o cabeçalho que o IEDB escreve."""
    rows = cc.conservancy_rows([("A_E_M_1_8", "ACDEFGHI")], PROTEINS, threshold, operator)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(cc.rows_to_csv(rows, cc.csv_columns(operator, threshold)), encoding="utf-8")
    return path


def poa1_dir(root, species=("A",)):
    """O diretório `Conservancy Analysis` do POA1, com um FASTA de epítopos por espécie."""
    root.mkdir(parents=True, exist_ok=True)
    for specie in species:
        (root / f"{specie}_epitopes.fasta").write_text(">A_E_M_1_8\nACDEFGHI\n", encoding="utf-8")
    return root


# ------------------------------------------------------------------ inventário do lote
def test_inventory_records_only_the_csvs_the_batch_produced(tmp_path):
    """
    Trava: write_metadata gravava um SHA-256 para TODO .csv do diretório, não só os do lote.

    O arquivo de anotações do pesquisador entrava no inventário; editá-lo passava a derrubar
    inspect_directory com "is missing or changed" — e com ele as etapas 4, 5 e o estágio do
    POA2 —, sem recuperação pela interface a não ser refazer a conservância inteira.
    """
    conservancy_csv(tmp_path / "A_conservancy.csv")
    notas = tmp_path / "notas.csv"
    notas.write_text("epitopo,comentario\nACDEFGHI,revisar\n", encoding="utf-8")

    cp.write_metadata(str(tmp_path), 70, ">=", source="local",
                      species={"A": {"csv": "A_conservancy.csv", "n_proteins": 2}})
    assert sorted(cp.read_metadata(str(tmp_path)).files) == ["A_conservancy.csv"]

    notas.write_text("epitopo,comentario\nACDEFGHI,descartar\n", encoding="utf-8")
    assert cp.inspect_directory(str(tmp_path)).origin == "meta"
    # O POA2 concatena todo CSV do diretório, então o arquivo estranho continua sendo reportado
    # — por um caminho que o usuário consegue resolver, em vez de um erro sem saída.
    assert cp.stale_csvs(str(tmp_path)) == ["notas.csv"]


def test_manual_declaration_leaves_a_csv_that_is_not_a_result_out_of_the_inventory(tmp_path):
    """
    Trava: o mesmo inventário indiscriminado na etapa 4 (CSVs enviados à mão).

    Ali não há lista de espécies e o diretório é o lote — mas um .csv sem a coluna de percentual
    não é um resultado de conservância, e o que foi produzido de verdade continua conferido.
    """
    alvo = conservancy_csv(tmp_path / "A_conservancy.csv")
    notas = tmp_path / "notas.csv"
    notas.write_text("epitopo,comentario\nACDEFGHI,revisar\n", encoding="utf-8")

    cp.write_metadata(str(tmp_path), 70, ">=", source="manual")
    assert sorted(cp.read_metadata(str(tmp_path)).files) == ["A_conservancy.csv"]

    notas.write_text("epitopo,comentario\nACDEFGHI,descartar\n", encoding="utf-8")
    assert cp.inspect_directory(str(tmp_path)).origin == "meta"

    alvo.write_text(alvo.read_text(encoding="utf-8").replace("100.00%", "50.00%"),
                    encoding="utf-8")
    with pytest.raises(ValueError, match="changed"):
        cp.inspect_directory(str(tmp_path))


# ------------------------------------------------------------------ lote sem resultado
def test_a_batch_that_produced_no_csv_declares_no_threshold(tmp_path):
    """
    Trava: a guarda `if written:` sumiu, e um lote sem nenhum CSV gravava conservancy_meta.json.

    O diretório passava a se declarar descrito (origin 'meta') com inventário vazio: a etapa 5
    deixava de pedir o limiar e o POA2 escrevia uma planilha vazia rotulada com um limiar que
    não foi aplicado a dado nenhum. Sem metadado o diretório é 'unknown' — que é o que ele é —
    e a interface pergunta.
    """
    vazio = tmp_path / "poa1"
    vazio.mkdir()                                   # o POA1 filtrou todos os epítopos
    saida = tmp_path / "csvs"

    assert cc.run_conservancy_for_dir(str(vazio), PROTEINS, 70, str(saida)) == []
    assert not (saida / cp.METADATA_FILENAME).exists()

    fonte = cp.inspect_directory(str(saida))
    assert (fonte.origin, fonte.csv_count, fonte.threshold) == ("unknown", 0, None)


def test_an_empty_batch_does_not_inherit_the_previous_batch_threshold(tmp_path):
    """
    Trava: o mesmo defeito quando já havia um lote no diretório.

    A execução vazia arquivava os CSVs anteriores e instalava só o metadado, deixando o limiar
    declarado sem nenhum CSV que o sustentasse.
    """
    origem = poa1_dir(tmp_path / "poa1")
    saida = tmp_path / "csvs"
    assert cc.run_conservancy_for_dir(str(origem), PROTEINS, 70, str(saida))
    assert cp.inspect_directory(str(saida)).origin == "meta"

    (origem / "A_epitopes.fasta").unlink()          # nenhuma espécie sobrou
    assert cc.run_conservancy_for_dir(str(origem), PROTEINS, 70, str(saida)) == []
    assert cp.map_EpConservFiles(str(saida)) == []
    assert cp.inspect_directory(str(saida)).origin == "unknown"


# ------------------------------------------------------------------ operador por família
def test_header_and_metadata_operators_are_compared_by_family(tmp_path):
    """
    Trava: a checagem nova cruzava cabeçalho e metadado pelo SÍMBOLO do operador.

    O IEDB escreve '<=' onde o POA2 escreve '<'. Comparado ao pé da letra, um export legítimo
    derrubava inspect_directory — e com ele as etapas 4, 5, 7 e check_threshold — num diretório
    que o resto do módulo (operator_family, check_threshold) declara compatível. A comparação
    por família não pode, em troca, afrouxar o desacordo que a checagem existe para pegar.
    """
    iedb = tmp_path / "iedb"
    conservancy_csv(iedb / "A_conservancy.csv", 70, "<=")
    cp.write_metadata(str(iedb), 70, "<", source="iedb")

    assert cp.inspect_directory(str(iedb)).origin == "meta"
    assert cp.check_threshold(str(iedb), 70, "<") is None

    misturado = tmp_path / "misturado"
    conservancy_csv(misturado / "A_conservancy.csv", 70, ">=")
    cp.write_metadata(str(misturado), 70, "<", source="manual")
    with pytest.raises(ValueError, match="does not match"):
        cp.inspect_directory(str(misturado))


# ------------------------------------------------------------------ inventário ilegível
@pytest.mark.parametrize("inventario", [["A_conservancy.csv"], "A_conservancy.csv",
                                        {"A_conservancy.csv": 5}, 5])
def test_an_inventory_of_the_wrong_type_is_ignored_instead_of_crashing(tmp_path, inventario):
    """
    Trava: read_metadata copiava 'files' do JSON sem checar o tipo.

    Um conservancy_meta.json editado à mão (ou escrito por outra versão) escapava do próprio
    `except` e reaparecia como AttributeError em meta.files.items(), que os `except ValueError`
    da interface não pegam — a página do Streamlit mostrava traceback. E stale_csvs lia o
    mesmo valor como inventário, devolvendo uma resposta sem lastro.
    """
    conservancy_csv(tmp_path / "A_conservancy.csv")
    (tmp_path / cp.METADATA_FILENAME).write_text(json.dumps(
        {"threshold": 70, "operator": ">=", "source": "manual", "version": 2,
         "files": inventario}), encoding="utf-8")

    assert cp.read_metadata(str(tmp_path)) is None
    fonte = cp.inspect_directory(str(tmp_path))
    assert (fonte.origin, fonte.threshold, fonte.operator) == ("header", 70.0, ">=")
    assert cp.stale_csvs(str(tmp_path)) == []
