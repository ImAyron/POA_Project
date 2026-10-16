"""O algoritmo de conservância, fixado contra valores calculados à mão.

O IEDB não expõe o Epitope Conservancy Analysis por API (verificado em 12/10/2026), então a
equivalência numérica com o site não pode ser verificada automaticamente — ela depende de uma
exportação manual, comparada por ``tools/validate_conservancy.py``.

O que *pode* ser verificado sem o site, e é o que está aqui, é a fidelidade ao método publicado
(Bui et al., BMC Bioinformatics 2007): janela sem gaps do tamanho do epítopo, deslizada por toda
a proteína, identidade = casamentos/L, e o **máximo** sobre os offsets. Cada número abaixo foi
calculado à mão e está mostrado na forma de conta, não como constante mágica.

Estes testes travam o algoritmo. Se um deles quebrar, o cálculo mudou — e qualquer resultado
publicado com a versão anterior deixa de valer.
"""
import pytest

from poa.services import conservancy_client as cc


# --------------------------------------------------------------------------- a janela deslizante
def test_identity_is_matches_over_epitope_length():
    """
    G G A Y L G G  contra  A Y I, melhor janela no offset 2:

        A Y I
        A Y L   -> A==A, Y==Y, L!=I  ->  2/3
    """
    assert cc.max_identity("AYI", "GGAYLGG") == pytest.approx(2 / 3 * 100)


def test_the_maximum_wins_not_the_first_window():
    """
    'AYXAYI' tem duas janelas que casam parcialmente; vale a melhor, não a primeira.

        offset 0:  A Y X  vs  A Y I  ->  2/3 = 66.67%
        offset 3:  A Y I  vs  A Y I  ->  3/3 = 100%
    """
    assert cc.max_identity("AYI", "AYXAYI") == 100.0
    # Sem a segunda janela, o mesmo prefixo dá 2/3 — prova que o 100% veio do offset 3.
    assert cc.max_identity("AYI", "AYXAYQ") == pytest.approx(2 / 3 * 100)


def test_no_gaps_are_allowed():
    """
    Uma inserção na proteína não é absorvida: o método é ungapped.

        epitopo  A Y I
        proteina A Y Q I   -> nenhuma janela de 3 contem 'AYI'
        janelas: AYQ (2/3), YQI (1/3)  -> maximo 2/3
    """
    assert cc.max_identity("AYI", "AYQI") == pytest.approx(2 / 3 * 100)


def test_comparison_is_case_insensitive():
    assert cc.max_identity("ayi", "GGAYIGG") == 100.0
    assert cc.max_identity("AYI", "ggayigg") == 100.0


def test_epitope_as_long_as_the_protein_has_exactly_one_window():
    assert cc.max_identity("AYI", "AYI") == 100.0
    assert cc.max_identity("AYI", "AYQ") == pytest.approx(2 / 3 * 100)


# --------------------------------------------------------------------------- a fronteira do limiar
@pytest.mark.parametrize("operator, conta_como_match", [(">=", True), ("<", False)])
def test_identity_exactly_at_the_threshold(operator, conta_como_match):
    """
    O caso de fronteira, que decide a classificação de um epítopo e é fácil de errar.

    Epítopo de 4 resíduos com 3 casamentos = 75.0%, limiar 75:
      '>=' (conservado) inclui;  '<' (único) exclui.
    """
    assert cc.matches_threshold(75.0, 75, operator) is conta_como_match


def test_the_two_operators_partition_every_identity():
    """Nenhuma identidade pode contar para os dois objetivos, nem para nenhum."""
    for identity in (0.0, 33.3, 69.99, 70.0, 70.01, 100.0):
        ge = cc.matches_threshold(identity, 70, ">=")
        lt = cc.matches_threshold(identity, 70, "<")
        assert ge != lt, f"identidade {identity} classificada como {ge}/{lt}"


def test_unknown_operator_is_refused():
    with pytest.raises(ValueError):
        cc.matches_threshold(80.0, 70, "~=")


# --------------------------------------------------------------------------- a linha do CSV
def _linha(epitopo, proteinas, limiar, operador=">="):
    rows = cc.conservancy_rows([("SP_PROT_M_1_9", epitopo)], proteinas, limiar, operador)
    return rows[0]


def test_percent_column_is_matched_over_comparable():
    """
    Três proteínas, identidades 100%, 2/3 e 1/3 contra 'AYI'; limiar 70 conta só a primeira.
    """
    proteinas = [("P1", "GGAYIGG"), ("P2", "GGAYLGG"), ("P3", "GGAQLGG")]
    row = _linha("AYI", proteinas, 70)

    percent = row[cc.percent_column(">=", 70)]
    assert percent == "33.33% (1/3)"
    assert row["Maximum identity"] == "100.00%"           # a P1
    assert row["Minimum identity"] == "33.33%"            # a P3: so o G final casa


def test_minimum_and_maximum_are_over_proteins_not_over_windows():
    """
    'Minimum identity' é a proteína mais divergente do conjunto, não a pior janela de uma
    proteína — a pior janela é quase sempre 0% e a coluna não teria informação nenhuma.
    """
    proteinas = [("P1", "AYI"), ("P2", "AYQ")]
    row = _linha("AYI", proteinas, 70)

    assert row["Maximum identity"] == "100.00%"
    assert row["Minimum identity"] == f"{2 / 3 * 100:.2f}%"


def test_unique_objective_counts_the_proteins_below_the_threshold():
    """Com '<' o papel se inverte: conta quem NÃO contém o epítopo."""
    proteinas = [("P1", "GGAYIGG"), ("P2", "GGAQLGG")]
    row = _linha("AYI", proteinas, 70, operador="<")

    # P1 = 100% (nao conta), P2 = 33.33% (conta)
    assert row[cc.percent_column("<", 70)] == "50.00% (1/2)"


def test_proteins_too_short_leave_the_denominator():
    """
    O conjunto tem 3 proteínas, mas só 2 podem conter um epítopo de 3 resíduos.

    A terceira não é "0% idêntica": ela não é comparável. Contá-la no denominador diluiria a
    conservância com não-evidência, e sob '<' a contaria como prova de unicidade.
    """
    proteinas = [("P1", "GGAYIGG"), ("P2", "GGAYLGG"), ("CURTA", "GG")]

    conservado = _linha("AYI", proteinas, 70)
    assert conservado[cc.percent_column(">=", 70)] == "50.00% (1/2)"
    assert conservado["Minimum identity"] == f"{2 / 3 * 100:.2f}%"   # a P2, nao 0.00%

    unico = _linha("AYI", proteinas, 70, operador="<")
    assert unico[cc.percent_column("<", 70)] == "50.00% (1/2)"       # a P2 so; a CURTA nao conta


def test_epitope_length_column_matches_the_sequence():
    row = _linha("AYIQ", [("P1", "GGAYIQGG")], 70)
    assert row["Epitope length"] == 4


def test_single_protein_can_only_give_zero_or_one_hundred():
    """
    A razão de o aviso de conjunto unitário existir: a coluna é binária e nenhum dos dois
    valores significa conservação.
    """
    identica = _linha("AYI", [("P1", "GGAYIGG")], 70)
    distante = _linha("AYI", [("P1", "GGAQLGG")], 70)

    assert identica[cc.percent_column(">=", 70)] == "100.00% (1/1)"
    assert distante[cc.percent_column(">=", 70)] == "0.00% (0/1)"


def test_empty_comparison_set_does_not_raise():
    row = _linha("AYI", [], 70)
    assert row[cc.percent_column(">=", 70)] == "0.00% (0/0)"
