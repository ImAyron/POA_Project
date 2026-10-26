"""Proteína e espécie alinhadas ao `-f`, registro a registro.

Dois casos reais, com exigências opostas, no mesmo pipeline:

* um consenso de E2 de alphavirus com cabeçalho ``alphavirus_e2_consenso`` — espécie *alphavirus*,
  proteína *E2* — e o JSON do BepiPred-2.0 com a mesma string. Enviados pelo usuário, trazem a
  ordem da ferramenta que os gerou e **precisam** ser invertidos;
* o `.txt` que a própria interface gera pelo EMBOSS antigenic, produzido **a partir do `-f` já
  normalizado**, que sai na ordem da convenção e **não pode** ser invertido.

Uma inversão em bloco corrigia o primeiro e estragava o segundo. Comparar cada epítopo com o `-f`
resolve os dois sem ninguém declarar nada: a etapa 1 estabelece a verdade, e o POA1 alinha o que
discordar dela.
"""
import json

from poa.core import pipeline
from poa.gui import backend


def _b2_json(path, key="alphavirus_e2_consenso"):
    """Um JSON do BepiPred-2.0 com uma chave de antígeno e dois resíduos antigênicos."""
    path.write_text(json.dumps({"antigens": {key: {
        "AA": list("MKTAYIA"), "PRED": [0, 0.9, 0.9, 0, 0, 0, 0]}}}), encoding="utf-8")
    return str(path)


def _pap_txt(path, header="E2_ALPHAVIRUS_CONSENSO"):
    """O que o EMBOSS antigenic produz: cabeçalho vindo do próprio -f, já na ordem certa."""
    path.write_text(f">{header}\nn\tStart\tSequence\tEnd\n1\t2\tKTAYI\t6\n", encoding="utf-8")
    return str(path)


def _reference(tmp_path, header="E2_ALPHAVIRUS_CONSENSO"):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">{header}\nMKTAYIA\n", encoding="utf-8")
    return str(fasta)


def _run(tmp_path, files, reference, name="out"):
    args = backend.build_poa1_args(files, {}, reference, str(tmp_path / name))
    return pipeline.run_poa1(args)


# --------------------------------------------------------------- o alinhamento
def test_an_uploaded_file_in_the_other_order_is_aligned(tmp_path):
    """
    O caso original: o JSON do usuário traz `alphavirus_e2`, o `-f` diz E2/ALPHAVIRUS.

    Sem isso, o `-f` e os epítopos da mesma string discordavam, e a etapa 3 relatava zero
    coincidência entre eles.
    """
    result = _run(tmp_path, {"b2": _b2_json(tmp_path / "ep.json")}, _reference(tmp_path))

    assert {s.upper() for s in result.predictions["Specie"]} == {"ALPHAVIRUS"}
    assert {s.upper() for s in result.predictions["Protein"]} == {"E2"}


def test_a_file_the_app_generated_is_left_alone(tmp_path):
    """
    O caso que a inversão em bloco quebrava.

    O EMBOSS antigenic roda sobre `ss.proteins_path` — o `-f` **já normalizado** — então o `.txt`
    sai com `E2_ALPHAVIRUS_CONSENSO`. Inverter por decreto o tornava errado; comparar com o `-f`
    o reconhece como já correto.
    """
    result = _run(tmp_path, {"p": _pap_txt(tmp_path / "pap.txt")}, _reference(tmp_path))

    assert {s.upper() for s in result.predictions["Specie"]} == {"ALPHAVIRUS"}
    assert {s.upper() for s in result.predictions["Protein"]} == {"E2"}


def test_the_two_kinds_coexist_in_one_run(tmp_path):
    """
    O que torna a detecção necessária em vez de conveniente.

    Numa mesma execução há arquivos enviados e arquivos gerados, com ordens opostas. Nenhuma
    decisão única — inverter tudo ou não inverter nada — acerta os dois.
    """
    result = _run(tmp_path, {"b2": _b2_json(tmp_path / "ep.json"),
                             "p": _pap_txt(tmp_path / "pap.txt")}, _reference(tmp_path))

    assert len(result.predictions) > 1, "os dois métodos precisam ter contribuído"
    assert {s.upper() for s in result.predictions["Specie"]} == {"ALPHAVIRUS"}
    assert {s.upper() for s in result.predictions["Protein"]} == {"E2"}


def test_a_pair_the_reference_does_not_have_is_left_untouched(tmp_path):
    """
    Espécie ausente do `-f` já é relatada pelo diagnóstico da etapa 3.

    Adivinhar aqui só faria a resposta errada parecer deliberada.
    """
    b2 = _b2_json(tmp_path / "ep.json", key="ns1_denv1_ref")
    result = _run(tmp_path, {"b2": b2}, _reference(tmp_path))

    assert {s.upper() for s in result.predictions["Specie"]} == {"DENV1"}
    assert {s.upper() for s in result.predictions["Protein"]} == {"NS1"}


def test_the_per_species_fasta_uses_the_aligned_species(tmp_path):
    """Os arquivos da etapa 4 são nomeados pela espécie, e é por eles que a conservância agrupa."""
    from pathlib import Path

    result = _run(tmp_path, {"b2": _b2_json(tmp_path / "ep.json")}, _reference(tmp_path))
    written = sorted(p.name for p in Path(result.conservancy_dir).glob("*"))

    assert any(n.upper().startswith("ALPHAVIRUS") for n in written), written
    assert not any(n.upper().startswith("E2_") for n in written), written


def test_the_reference_and_the_predictions_end_up_agreeing(tmp_path):
    """O sintoma que o usuário via: zero coincidência entre o -f e os epítopos."""
    reference = _reference(tmp_path)
    result = _run(tmp_path, {"b2": _b2_json(tmp_path / "ep.json")}, reference)

    diag = backend.species_diagnostics(reference, result.predictions)

    assert diag["predicted"] == ["ALPHAVIRUS"]
    assert not diag["missing_from_reference"], diag


# --------------------------------------------------------------- as peças, isoladas
def test_reference_pairs_are_read_uppercased(tmp_path):
    """Os parsers diferem na caixa: o b2 preserva a do JSON, um FASTA sobe para maiúscula."""
    reference = _reference(tmp_path, header="e2_alphavirus_consenso")

    assert pipeline.reference_header_pairs(reference) == {("E2", "ALPHAVIRUS")}


def test_a_reference_without_a_second_field_yields_no_pairs(tmp_path):
    """Sem segundo campo não há par, e sem par nada é alinhado — nem para o lado errado."""
    reference = _reference(tmp_path, header="semseparador")

    assert pipeline.reference_header_pairs(reference) == set()


def test_align_reports_how_many_it_moved(tmp_path):
    """A contagem vai para o log: uma correção silenciosa é indistinguível de nada ter havido."""
    import pandas as pd

    from poa.core.models import STANDARD_COLUMNS

    frame = pd.DataFrame([
        {"Method": "m", "Specie": "E2", "Protein": "ALPHAVIRUS", "ID_Sequence": "x",
         "Initial Position": 1, "Final Position": 3, "Peptide Sequence": "MKT"},
        {"Method": "m", "Specie": "ALPHAVIRUS", "Protein": "E2", "ID_Sequence": "x",
         "Initial Position": 4, "Final Position": 6, "Peptide Sequence": "AYI"},
        {"Method": "m", "Specie": "ZIKA", "Protein": "NS1", "ID_Sequence": "x",
         "Initial Position": 7, "Final Position": 9, "Peptide Sequence": "AKQ"},
    ])[list(STANDARD_COLUMNS)]

    out, aligned, unresolved = pipeline.align_to_reference(frame, _reference(tmp_path))

    assert aligned == 1 and unresolved == 1
    assert list(out["Specie"]) == ["ALPHAVIRUS", "ALPHAVIRUS", "ZIKA"]
    assert list(out.columns) == list(STANDARD_COLUMNS)


def test_an_empty_table_is_returned_unchanged(tmp_path):
    """Nenhum epítopo sobrevivendo aos filtros é caso normal, não erro."""
    import pandas as pd

    out, aligned, unresolved = pipeline.align_to_reference(pd.DataFrame(), _reference(tmp_path))

    assert out.empty and aligned == 0 and unresolved == 0


# --------------------------------------------------------------- a prévia da etapa 2
def test_the_step2_preview_shows_the_order_poa1_will_use(tmp_path):
    """
    A prévia lê o arquivo direto, sem passar pelo POA1, e precisava da mesma regra.

    Era o sintoma que restava: a etapa 1 concordava com o usuário e a etapa 2 logo adiante ainda
    anunciava `Espécies: E2` para uma proteína chamada E2.
    """
    b2 = _b2_json(tmp_path / "ep.json")
    reference = _reference(tmp_path)

    assert backend.summarize_sources({"b2": [b2]})[0].species == ["E2"]
    assert backend.summarize_sources({"b2": [b2]}, reference=reference)[0].species == ["ALPHAVIRUS"]


def test_the_preview_leaves_a_generated_file_alone(tmp_path):
    """E não pode inverter o que já está certo — o mesmo erro da inversão em bloco."""
    fasta = tmp_path / "b3.fasta"
    fasta.write_text(">E2_ALPHAVIRUS_CONSENSO\nmktAYIAkqrq\n", encoding="utf-8")

    preview = backend.summarize_sources({"b3": [str(fasta)]}, reference=_reference(tmp_path))[0]

    assert preview.species == ["ALPHAVIRUS"] and preview.proteins == ["E2"]
    assert preview.rows[0]["Espécie"] == "ALPHAVIRUS"
    assert preview.rows[0]["Proteína"] == "E2"


def test_one_species_by_two_methods_is_one_file(tmp_path):
    """
    Os parsers discordam na caixa, e a etapa 4 nomeia os arquivos por esta coluna.

    O b2 preserva o `alphavirus` do JSON e o PAP/IMED sobe para `ALPHAVIRUS`: a mesma espécie
    chegando por dois métodos era escrita como dois arquivos, e a conservância comparava cada
    metade contra um conjunto como se fossem organismos diferentes. Um registro reconhecido no
    `-f` passa a ser reescrito com a grafia dele, não só com a ordem.
    """
    from pathlib import Path

    result = _run(tmp_path, {"b2": _b2_json(tmp_path / "ep.json"),
                             "p": _pap_txt(tmp_path / "pap.txt")}, _reference(tmp_path))
    written = sorted(p.name for p in Path(result.conservancy_dir).glob("*_epitopes.fasta"))

    assert written == ["ALPHAVIRUS_epitopes.fasta"], written
    assert len({str(s) for s in result.predictions["Specie"]}) == 1
