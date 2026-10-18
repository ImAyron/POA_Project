"""Reabrir uma análise devolve cada arquivo à etapa em que foi usado.

Antes, reabrir restaurava só os CSVs de conservância e a referência, atrás de um botão na etapa 5.
Todo o resto — os arquivos de predição que as pessoas enviam à mão na etapa 2, os conjuntos de
comparação por espécie da etapa 4, a estrutura da etapa 7 — ficava no disco enquanto a interface
mostrava campos vazios, então continuar uma análise significava reenviar o que ela já continha.

O nome do arquivo não diz de que etapa ele veio: um FASTA de referência e um de "outros
preditores" são ambos FASTAs em ``inputs/``, e um arquivo que a pessoa nomeou não segue convenção
nenhuma. Por isso a procedência é **gravada na chegada**; a classificação por nome existe só para
as análises que já estavam no disco antes desse registro.
"""
import json
import logging
import os
from pathlib import Path

import pandas as pd
import pytest

from poa.gui import backend


@pytest.fixture
def ctx(tmp_path, monkeypatch):
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    return backend.create_run(tmp_path, "teste restore")


def _upload(ctx, dest, name, text, **provenance):
    return backend.save_upload(dest, name, text.encode("utf-8"), ctx=ctx, **provenance)


# --------------------------------------------------------------- o histórico
def test_history_records_where_each_file_came_from(ctx):
    _upload(ctx, ctx.inputs_dir, "ref_0_meu.fasta", ">NS1_DENV1_ref\nMKT\n",
            role="reference", original_name="meu.fasta")
    _upload(ctx, ctx.inputs_dir, "others_0_x.fasta", ">E_DENV1_OTHER_1_3\nMKT\n",
            role="prediction", method="x", original_name="x.fasta")
    _upload(ctx, ctx.inputs_dir, "comparison_DENV1.fasta", ">E_DENV1_a\nMKT\n",
            role="comparison_set", specie="DENV1", original_name="world.fasta")

    history = backend.read_input_history(ctx)
    by_path = {item["path"]: item for item in history}

    assert by_path["inputs/ref_0_meu.fasta"]["role"] == "reference"
    assert by_path["inputs/others_0_x.fasta"]["method"] == "x"
    assert by_path["inputs/comparison_DENV1.fasta"]["specie"] == "DENV1"
    # o nome que a pessoa deu é guardado, não só o nome interno
    assert by_path["inputs/comparison_DENV1.fasta"]["original_name"] == "world.fasta"


def test_history_paths_are_relative_to_the_analysis_folder(ctx):
    """
    A mesma pasta é aberta do Windows e da WSL, onde o caminho absoluto é outro.

    Gravar o caminho absoluto de um lado tornaria o registro inútil do outro — e esse vai-e-vem
    entre os dois ambientes é justamente o fluxo de trabalho do projeto.
    """
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">X_SP_r\nMKT\n", role="reference")

    recorded = backend.read_input_history(ctx)[0]["path"]

    assert recorded == "inputs/ref_0_a.fasta"
    assert ":" not in recorded and not recorded.startswith("/")


def test_re_uploading_the_same_path_does_not_duplicate_the_entry(ctx):
    for _ in range(3):
        _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">X_SP_r\nMKT\n", role="reference")

    assert len(backend.read_input_history(ctx)) == 1


def test_an_unreadable_history_is_ignored_not_fatal(ctx):
    (ctx.root / backend.INPUT_HISTORY_FILENAME).write_text("{nao é json", encoding="utf-8")

    assert backend.read_input_history(ctx) == []
    # e gravar de novo recompõe o arquivo em vez de estourar
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">X_SP_r\nMKT\n", role="reference")
    assert len(backend.read_input_history(ctx)) == 1


def test_saving_without_a_role_writes_the_file_and_no_history(ctx):
    """O histórico é um acréscimo, nunca uma condição para salvar."""
    path = backend.save_upload(ctx.inputs_dir, "solto.txt", b"dados")

    assert backend.read_input_history(ctx) == []
    assert open(path, encoding="utf-8").read() == "dados"


# --------------------------------------------------------------- a restauração
def test_every_file_goes_back_to_its_own_step(ctx):
    _upload(ctx, ctx.inputs_dir, "ref_0_prot.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")
    _upload(ctx, ctx.inputs_dir, "bepipred3_0_b.fasta", ">NS1_DENV1_ref\nmktAYIa\n",
            role="prediction", method="b3")
    _upload(ctx, ctx.inputs_dir, "netctl_0_n.html", "<html></html>",
            role="prediction", method="n")
    _upload(ctx, ctx.inputs_dir, "comparison_DENV1.fasta", ">NS1_DENV1_w\nMKTAYIA\n",
            role="comparison_set", specie="DENV1")
    _upload(ctx, ctx.inputs_dir, "structure_4C23.pdb", "ATOM\nEND\n", role="structure")
    _upload(ctx, ctx.mhcii_dir, "NS1_DENV1.html", "<html></html>",
            role="prediction", method="m")

    restored = backend.restore_analysis(ctx)

    assert restored.source == "history"
    assert [p.endswith("ref_0_prot.fasta") for p in restored.references] == [True]
    assert restored.prediction_files["b3"][0].endswith("bepipred3_0_b.fasta")
    assert restored.prediction_files["n"][0].endswith("netctl_0_n.html")
    assert restored.prediction_files["m"] == [str(ctx.mhcii_dir)]   # MHC-II é a pasta, não o arquivo
    assert restored.comparison_sets["DENV1"].endswith("comparison_DENV1.fasta")
    assert restored.structure.endswith("structure_4C23.pdb")
    assert restored.anything


def test_a_file_deleted_from_disk_is_dropped_from_the_restore(ctx):
    """O histórico registra o passado; o disco decide o presente."""
    path = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">X_SP_r\nMKT\n", role="reference")
    (ctx.inputs_dir / "ref_0_a.fasta").unlink()

    restored = backend.restore_analysis(ctx)

    assert restored.references == []
    assert backend.read_input_history(ctx), "a linha continua no histórico"
    assert path  # o caminho existia quando foi gravado


def test_an_empty_analysis_restores_nothing(ctx):
    restored = backend.restore_analysis(ctx)

    assert not restored.anything
    assert restored.summary() == []


# ------------------------------------------- as duas contas de procedência, lidas juntas
def test_a_new_upload_does_not_hide_what_the_analysis_already_had(ctx):
    """
    O registro e o nome cobrem conjuntos diferentes de arquivos, então são lidos juntos.

    Uma análise anterior ao registro ganha histórico assim que recebe um upload qualquer. Tratar
    as duas contas como exclusivas fazia uma única linha desligar a varredura inteira: tudo o que
    já estava no disco saía do campo de visão da reabertura, a pessoa reenviava o que já estava
    lá, e o aviso de "classificado pelo nome" sumia justo quando a restauração passou a fazer
    menos.
    """
    (ctx.inputs_dir / "ref_0_prot.fasta").write_text(">NS1_DENV1_ref\nMKTAYIA\n", encoding="utf-8")
    (ctx.inputs_dir / "others_0_x.fasta").write_text(">NS1_DENV1_O_1_3\nMKT\n", encoding="utf-8")
    (ctx.inputs_dir / "comparison_DENV1.fasta").write_text(">NS1_DENV1_w\nMKT\n", encoding="utf-8")
    (ctx.inputs_dir / "structure_1ABC.pdb").write_text("ATOM\nEND\n", encoding="utf-8")
    _upload(ctx, ctx.inputs_dir, "comparison_DENV2.fasta", ">NS1_DENV2_w\nMKT\n",
            role="comparison_set", specie="DENV2")

    restored = backend.restore_analysis(ctx)

    assert restored.references[0].endswith("ref_0_prot.fasta")
    assert restored.prediction_files["x"][0].endswith("others_0_x.fasta")
    assert restored.structure.endswith("structure_1ABC.pdb")
    assert sorted(restored.comparison_sets) == ["DENV1", "DENV2"]
    assert restored.source == "mixed"


def test_the_record_wins_the_single_slots_over_a_name(ctx):
    """
    Onde as duas contas preenchem o mesmo campo único, o que foi gravado vence o palpite.

    A varredura é lida primeiro justamente por isso: a estrutura é um campo só, e uma estrutura
    antiga no disco não pode sobrescrever a que a pessoa enviou — inverter a ordem reintroduz o
    erro sem quebrar mais nada.
    """
    (ctx.inputs_dir / "structure_ANTIGA.pdb").write_text("ATOM\nEND\n", encoding="utf-8")
    _upload(ctx, ctx.inputs_dir, "structure_NOVA.pdb", "ATOM\nEND\n", role="structure")

    restored = backend.restore_analysis(ctx)

    assert restored.structure.endswith("structure_NOVA.pdb")
    assert restored.source == "mixed"


def test_a_merge_assembled_for_poa1_is_not_restored_beside_its_own_parts(ctx):
    """
    O arquivo que ``consolidate_sources`` monta não é entrada de ninguém.

    Ele carrega o mesmo radical das partes ('bepipred3_merged.fasta' ao lado de
    'bepipred3_0_DENV1.fasta'), então uma varredura ingênua o ofereceria ao POA1 junto com as
    partes de que ele mesmo foi feito, e cada epítopo seria contado duas vezes.
    """
    for i, sp in enumerate(["DENV1", "DENV2"]):
        _upload(ctx, ctx.inputs_dir, f"bepipred3_{i}_{sp}.fasta", f">NS1_{sp}_ref\nMKTAYIA\n",
                role="prediction", method="b3")
    parts = [str(p) for p in sorted(ctx.inputs_dir.glob("bepipred3_?_*.fasta"))]
    backend.consolidate_sources({"b3": parts}, ctx.inputs_dir)
    assert (ctx.inputs_dir / "bepipred3_merged.fasta").is_file()

    restored = backend.restore_analysis(ctx)

    assert len(restored.prediction_files["b3"]) == 2
    assert restored.source == "history"


def test_a_file_the_name_cannot_place_does_not_raise_the_warning(ctx):
    """
    Um arquivo que o nome não coloca em etapa nenhuma não é palpite sobre nada.

    Ele não restaura coisa alguma, então não há o que conferir: contá-lo faria uma análise de
    registro completo pedir conferência à toa.
    """
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">X_SP_r\nMKT\n", role="reference")
    (ctx.inputs_dir / "POA0_relatorio.json").write_text("{}", encoding="utf-8")

    assert backend.restore_analysis(ctx).source == "history"


# --------------------------------------------------------------- análises anteriores ao registro
@pytest.mark.parametrize("name, role, method", [
    ("ref_0_meu.fasta", "reference", ""),
    ("proteins_all.fasta", "merged_reference", ""),
    ("bepipred3_0_x.fasta", "prediction", "b3"),
    ("pap_antigenic.txt", "prediction", "p"),
    ("netctl_0_x.html", "prediction", "n"),
    ("others_0_x.fasta", "prediction", "x"),
    ("comparison_DENV2.fasta", "comparison_set", ""),
    ("structure_1ABC.pdb", "structure", ""),
    ("raw_dv1_dv2.json", "bepipred2_raw", ""),
    ("dv1_b2.json", "bepipred2_adapted", "b2"),
    ("dv1_ref.fasta", "bepipred2_reference", ""),
    ("E2_OROV.fasta", "reference", ""),          # saída do POA0, sem prefixo próprio
])
def test_legacy_names_are_classified(name, role, method):
    assert backend._classify_legacy_input(name) == (role, method)


@pytest.mark.parametrize("name, role, method", [
    ("ref_0_E_DENV1_ref.fasta", "reference", ""),   # o ID que a ajuda da etapa 1 sugere é 'ref'
    ("others_0_x_ref.fasta", "prediction", "x"),
    ("bepipred3_0_saida_ref.fasta", "prediction", "b3"),
    ("raw_denv1_b2.json", "bepipred2_raw", ""),
])
def test_the_prefix_app_py_wrote_beats_a_suffix_from_the_users_own_name(name, role, method):
    """
    O radical é a parte que o app escreve; tudo depois dele é o nome que a pessoa deu.

    Os sufixos '_b2.json' e '_ref.fasta' pertencem a um escritor só, o adaptador, cujo próprio
    radical é uma espécie digitada — a evidência mais fraca das duas. Testá-los primeiro deixava o
    nome do usuário sequestrar a classificação de qualquer arquivo que o app tivesse prefixado, e
    como nem 'bepipred2_reference' nem 'bepipred2_raw' são restaurados, o arquivo era descartado
    em silêncio com um resumo plausível e errado.
    """
    assert backend._classify_legacy_input(name) == (role, method)


@pytest.mark.parametrize("name, role, method", [
    ("PAPILLOMA_b2.json", "bepipred2_adapted", "b2"),
    ("PAPILLOMA_ref.fasta", "bepipred2_reference", ""),
    ("OTHERSVIRUS_b2.json", "bepipred2_adapted", "b2"),
])
def test_a_species_that_merely_starts_like_a_save_site_is_not_a_prediction(name, role, method):
    """
    Cada radical casa só até o separador que o app de fato escreve.

    'pap', 'others', 'netctl' e 'bepipred3' estão na tabela como radicais soltos: sem exigir o
    '_' ou o '.', inverter a ordem faria 'PAPILLOMA_b2.json' virar predição de antigenicidade.
    As duas mudanças são uma correção só.
    """
    assert backend._classify_legacy_input(name) == (role, method)


def test_a_legacy_reference_named_ref_fasta_comes_back_to_step_1(ctx):
    """A consequência do sufixo testado antes do radical, e não só a tabela: etapas perdidas."""
    (ctx.inputs_dir / "ref_0_E_DENV1_ref.fasta").write_text(">E_DENV1_ref\nMKTAYIA\n",
                                                            encoding="utf-8")
    (ctx.inputs_dir / "others_0_x_ref.fasta").write_text(">E_DENV1_O_1_3\nMKT\n", encoding="utf-8")

    restored = backend.restore_analysis(ctx)

    assert restored.references[0].endswith("ref_0_E_DENV1_ref.fasta")
    assert restored.prediction_files["x"][0].endswith("others_0_x_ref.fasta")
    assert restored.source == "inferred"


def test_the_first_guis_comparison_fasta_is_not_merged_into_the_reference(ctx):
    """
    'comparison.fasta' é o nome que a primeira GUI dava ao conjunto de comparação.

    Sem o separador no casamento ele caía no fallback de FASTA e entrava no -f da etapa 1 — um
    conjunto de diversidade somado às proteínas de referência, o pior dos dois erros. Sob a chave
    "COMPARISON" ele não casa com espécie nenhuma e fica inerte, aparecendo só no resumo.
    """
    (ctx.inputs_dir / "comparison.fasta").write_text(">NS1_DENV1_w\nMKT\n", encoding="utf-8")

    restored = backend.restore_analysis(ctx)

    assert restored.references == []
    assert [Path(p).name for p in restored.comparison_sets.values()] == ["comparison.fasta"]


def test_an_analysis_without_history_is_restored_by_name_and_says_so(ctx):
    """
    As análises que já estão no disco não têm registro de procedência.

    Elas são a razão de a classificação por nome existir: restaurá-las é o pedido, e um manifesto
    que só vale daqui pra frente não as alcançaria.
    """
    (ctx.inputs_dir / "ref_0_prot.fasta").write_text(">NS1_DENV1_ref\nMKTAYIA\n", encoding="utf-8")
    (ctx.inputs_dir / "others_0_x.fasta").write_text(">E_DENV1_O_1_3\nMKT\n", encoding="utf-8")
    assert not (ctx.root / backend.INPUT_HISTORY_FILENAME).exists()

    restored = backend.restore_analysis(ctx)

    assert restored.source == "inferred"
    assert restored.references[0].endswith("ref_0_prot.fasta")
    assert restored.prediction_files["x"][0].endswith("others_0_x.fasta")


def test_history_wins_over_the_name_when_both_exist(ctx):
    """
    Um arquivo que a pessoa nomeou não segue convenção nenhuma — o registro é que vale.

    Aqui o nome diria "predição de outros preditores"; o registro diz que entrou como referência.
    """
    _upload(ctx, ctx.inputs_dir, "others_0_x.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
            role="reference", original_name="minhas_proteinas.fasta")

    restored = backend.restore_analysis(ctx)

    assert restored.references[0].endswith("others_0_x.fasta")
    assert "x" not in restored.prediction_files


def test_summary_names_the_steps_the_files_belong_to(ctx):
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">X_SP_r\nMKT\n", role="reference")
    _upload(ctx, ctx.inputs_dir, "others_0_x.fasta", ">X_SP_O_1_3\nMKT\n",
            role="prediction", method="x")
    _upload(ctx, ctx.inputs_dir, "structure_1ABC.pdb", "ATOM\nEND\n", role="structure")

    lines = " | ".join(backend.restore_analysis(ctx).summary())

    assert "Etapa 1" in lines and "Etapa 2" in lines and "Etapa 7" in lines


# --------------------------------------------------------------- o BepiPred-2.0 da etapa 2
def _bepipred2_json(**antigens):
    """Um export do BepiPred-2.0, uma chave de antígeno por sequência."""
    return json.dumps({
        "info": {"failedjobs": 0, "size": len(antigens)},
        "antigens": {key: {"AA": list(seq), "BepiPred-2.0 score": [0.1] * len(seq)}
                     for key, seq in antigens.items()},
    }).encode("utf-8")


def _step_two_import(ctx, filename, data, **kwargs):
    """
    A etapa 2 inteira, como ``app.py`` a executa: grava o JSON cru e só então adapta.

    O envio tem de passar pelo ``save_upload`` antes — é ele que dá histórico à análise, e é com
    histórico que o defeito aparece. Um teste que chamasse ``import_bepipred2`` numa análise sem
    nenhum registro cairia na classificação por nome e passaria mesmo sem a correção.
    """
    backend.save_upload(ctx.inputs_dir, f"raw_{filename}", data, ctx=ctx,
                        role="bepipred2_raw", original_name=filename)
    return backend.import_bepipred2(ctx, filename, data, **kwargs)


def test_the_adapted_bepipred2_pair_comes_back(ctx):
    """
    O adaptador escreve seus dois arquivos sozinho, por fora do ``save_upload``.

    Ninguém os registrava, e o registro vence os nomes no disco; como o JSON cru é gravado logo
    antes de adaptar, toda importação real da etapa 2 já tem histórico. A etapa 2 reabria vazia,
    pedindo de novo um JSON que a própria análise já adaptou.
    """
    prep = _step_two_import(ctx, "bepipred_dv1.json",
                            _bepipred2_json(Sequence="MKTAYIA"), specie="DENV1")

    restored = backend.restore_analysis(ctx)

    assert restored.source == "history"
    assert len(restored.prepared_b2) == 1
    assert restored.prepared_b2[0].b2_json == prep.b2_json
    assert restored.prepared_b2[0].reference_fasta == prep.reference_fasta
    assert [e.header for e in restored.prepared_b2[0].entries] == ["E_DENV1_ref"]


def test_the_rebuilt_reference_is_not_restored_as_a_step_1_fasta(ctx):
    """
    O -f reconstruído volta dentro do par adaptado, nunca como um FASTA da etapa 1.

    Ele já é somado ao -f por ``prepared_b2[].reference_fasta``; devolvê-lo também como referência
    mesclaria a mesma proteína duas vezes em proteins_all.fasta. É por isso que o papel
    'bepipred2_reference' é gravado como procedência e não tem ramo na restauração.
    """
    _step_two_import(ctx, "b.json", _bepipred2_json(Sequence="MKTAYIA"), specie="DENV1")

    restored = backend.restore_analysis(ctx)

    assert restored.references == []
    assert restored.prepared_b2[0].reference_fasta.endswith("DENV1_ref.fasta")


def test_a_json_mapped_to_several_species_is_restored_as_one_record(ctx):
    """Uma submissão de várias espécies é um par só, pareado pelo radical 'DENV1-DENV2'."""
    _step_two_import(ctx, "b.json", _bepipred2_json(denv1="MKTAYIA", denv2="QRQMKTA"),
                     specie="DENV1",
                     mapping={"denv1": ("E", "DENV1", "ref"), "denv2": ("E", "DENV2", "ref")})

    restored = backend.restore_analysis(ctx)

    assert len(restored.prepared_b2) == 1
    assert restored.prepared_b2[0].species == ["DENV1", "DENV2"]
    assert restored.prepared_b2[0].b2_json.endswith("DENV1-DENV2_b2.json")


def test_each_imported_species_comes_back_as_its_own_record(ctx):
    """Importações sucessivas se acumulam na etapa 2, e todas têm de voltar."""
    _step_two_import(ctx, "dv1.json", _bepipred2_json(Sequence="MKTAYIA"), specie="DENV1")
    _step_two_import(ctx, "chikv.json", _bepipred2_json(Sequence="QRQMKTA"),
                     specie="CHIKV", protein="E1")

    restored = backend.restore_analysis(ctx)

    assert sorted(e.header for pr in restored.prepared_b2
                  for e in pr.entries) == ["E1_CHIKV_ref", "E_DENV1_ref"]


def test_a_half_deleted_adapter_pair_leaves_nothing_behind(ctx):
    """
    Meio par não vira meia restauração — e o -f órfão não escorrega para a etapa 1.

    Com o par gravado, a restauração passa a procurar os dois arquivos; apagar o JSON não pode
    derrubar a reabertura nem fazer o FASTA que sobrou entrar como proteína de referência, porque
    ele foi reconstruído das previsões e não é um envio da etapa 1.
    """
    prep = _step_two_import(ctx, "b.json", _bepipred2_json(Sequence="MKTAYIA"), specie="DENV1")
    Path(prep.b2_json).unlink()

    restored = backend.restore_analysis(ctx)

    assert restored.prepared_b2 == []
    assert restored.references == []


def test_a_conforming_bepipred2_json_comes_back_as_the_b2_input(ctx):
    """
    Quando as chaves já seguem 'Proteína_Espécie_ID', a etapa 2 usa o próprio envio como -b2.

    O papel gravado continuava dizendo "bruto" e a restauração não tinha ramo para ele: o arquivo
    caía fora de todos os ``elif`` e era descartado em silêncio, então a reabertura de uma análise
    cujo único arquivo era esse parecia uma análise vazia.
    """
    path = backend.save_upload(ctx.inputs_dir, "raw_ok.json",
                               _bepipred2_json(NS1_DENV1_ref="MKTAYIA"),
                               ctx=ctx, role="bepipred2_raw", original_name="ok.json")

    restored = backend.restore_analysis(ctx)

    assert restored.prediction_files["b2"] == [path]
    assert any("Etapa 2 (b2)" in line for line in restored.summary())


def test_an_unmapped_bepipred2_json_is_not_offered_to_poa1(ctx):
    """
    Devolver um JSON não mapeado seria pior que não devolver nada.

    O POA1 não resolve espécie nenhuma a partir dele e todos os epítopos ficam sob ``nan``, que é
    exatamente o que o adaptador existe para evitar — daí a conformidade ser decidida pelas
    chaves, com a mesma regra que o envio usa, e não pelo papel gravado.
    """
    backend.save_upload(ctx.inputs_dir, "raw_gen.json", _bepipred2_json(denv1="MKTAYIA"),
                        ctx=ctx, role="bepipred2_raw", original_name="gen.json")

    assert "b2" not in backend.restore_analysis(ctx).prediction_files


def test_a_legacy_raw_json_is_decided_by_its_antigen_keys(ctx):
    """
    O nome não diz se o JSON serve; o conteúdo diz, e o prefixo 'raw_' é anterior ao registro.

    Há análises assim no disco, e elas são justamente o público da classificação por nome.
    """
    (ctx.inputs_dir / "raw_ok.json").write_bytes(_bepipred2_json(NS1_DENV1_ref="MKTAYIA"))
    (ctx.inputs_dir / "raw_generico.json").write_bytes(_bepipred2_json(denv1="MKTAYIA"))

    restored = backend.restore_analysis(ctx)

    assert restored.source == "inferred"
    assert [Path(p).name for p in restored.prediction_files["b2"]] == ["raw_ok.json"]


def test_an_unreadable_bepipred2_json_does_not_break_the_restore(ctx):
    """
    Nenhum arquivo que a reabertura apenas esperava encontrar pode derrubá-la.

    Um JSON com lista no topo levanta AttributeError na leitura das chaves, e ``{"antigens":
    null}`` levanta TypeError — nenhum dos dois é OSError, ValueError ou KeyError, então um
    ``except`` estreito deixaria a análise inteira deixar de abrir por causa de um arquivo.
    """
    (ctx.inputs_dir / "raw_quebrado.json").write_text("{não é json", encoding="utf-8")
    (ctx.inputs_dir / "raw_sem_antigenos.json").write_text('{"antigens": null}', encoding="utf-8")
    (ctx.inputs_dir / "raw_lista.json").write_text('[{"antigens": {}}]', encoding="utf-8")
    (ctx.inputs_dir / "ref_0_prot.fasta").write_text(">NS1_DENV1_ref\nMKTAYIA\n", encoding="utf-8")

    restored = backend.restore_analysis(ctx)

    assert "b2" not in restored.prediction_files
    assert restored.references[0].endswith("ref_0_prot.fasta")


# ------------------------------------------------- a seleção atual de cada etapa (step_state)
# O histórico registra CHEGADAS: que um arquivo entrou, por qual etapa e com que nome. É um fato
# sobre o arquivo, e nunca deixa de ser verdade. Várias etapas, porém, guardam uma SELEÇÃO — o
# uploader de cada método SUBSTITUI a lista daquele método, BepiPred-2.0 e 3.0 se expulsam porque
# são o mesmo método de célula B, "Limpar espécies importadas" esvazia a etapa sem que arquivo
# nenhum chegue. Nenhuma dessas escolhas deixa rastro num registro de chegadas, e reabrir pelas
# chegadas ressuscitava o que a pessoa tinha abandonado. É a isso que step_state.json responde.
@pytest.fixture
def warned():
    """
    As linhas WARN do backend, colhidas no próprio logger.

    ``caplog`` não serve aqui: ``setup_logging`` desliga ``propagate`` no logger ``poa``, e o
    import de ``poa.gui.app`` o chama — então o handler que o pytest instala na raiz vê ou não vê
    as linhas conforme a ordem em que os arquivos de teste rodaram, e um aviso que só é conferido
    às vezes não confere nada.
    """
    messages = []

    class _Collect(logging.Handler):
        def emit(self, record):
            messages.append(record.getMessage())

    handler = _Collect(level=logging.WARNING)
    logger = logging.getLogger("poa.gui")
    saved_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.WARNING)
    yield messages
    logger.removeHandler(handler)
    logger.setLevel(saved_level)


def test_the_step_state_is_recorded_relative_to_the_analysis_folder(ctx):
    """
    Mesma razão do histórico: a pasta é aberta do Windows e da WSL, com caminhos diferentes.

    A seleção é lida de volta para dizer QUAIS arquivos a etapa tinha; gravada em absoluto de um
    lado, do outro ela aponta para caminhos que não existem e a etapa reabre vazia — o defeito que
    a camada inteira existe para impedir, com a aparência de um arquivo apagado.
    """
    kept = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
                   role="reference")
    backend.record_step_state(ctx, references=[kept], sources={"b3": [kept]}, ref_swap=True)

    state = json.loads((ctx.root / backend.STEP_STATE_FILENAME).read_text(encoding="utf-8"))

    assert state["references"] == ["inputs/ref_0_a.fasta"]
    # o aninhamento também: a seleção da etapa 2 é um dicionário de listas
    assert state["sources"] == {"b3": ["inputs/ref_0_a.fasta"]}
    assert state["ref_swap"] is True
    assert backend.read_step_state(ctx) == state


def test_recording_one_step_does_not_erase_what_another_recorded(ctx):
    """
    Cada etapa grava quando renderiza, e só o que é dela.

    As etapas 1, 2 e 5 escrevem no mesmo arquivo em momentos diferentes; se gravar fosse
    substituir, abrir a etapa 2 apagaria a seleção da etapa 1 que acabou de ser feita, e a pessoa
    perderia a escolha anterior ao andar para a frente.
    """
    backend.record_step_state(ctx, references=["inputs/a.fasta"], params_poa1={"bmin": 8})
    backend.record_step_state(ctx, sources={"b3": ["inputs/b.fasta"]})
    backend.record_step_state(ctx, params_poa2={"t": 90})

    state = backend.read_step_state(ctx)

    assert sorted(state) == ["params_poa1", "params_poa2", "references", "sources"]
    assert state["references"] == ["inputs/a.fasta"]


def test_recording_the_same_selection_again_does_not_rewrite_the_file(ctx):
    """
    Os pontos de gravação ficam em código que o Streamlit reexecuta a cada interação.

    A etapa 1 grava sua seleção toda vez que renderiza, isto é, a cada clique da análise inteira.
    Reescrever o arquivo quando nada mudou o reescreveria por nada, numa pasta cujo conteúdo é o
    registro do que a análise é.
    """
    backend.record_step_state(ctx, ref_swap=False, params_poa1={"bmin": 8})
    target = ctx.root / backend.STEP_STATE_FILENAME
    os.utime(target, (1_000_000, 1_000_000))

    backend.record_step_state(ctx, ref_swap=False, params_poa1={"bmin": 8})

    assert int(target.stat().st_mtime) == 1_000_000
    # e mudar de verdade volta a escrever
    backend.record_step_state(ctx, ref_swap=True)
    assert int(target.stat().st_mtime) != 1_000_000


@pytest.mark.parametrize("content", ["{não é json", "[1, 2, 3]", "null"])
def test_a_step_state_that_cannot_be_read_is_ignored_not_fatal(ctx, content):
    """
    Uma preferência ilegível vale menos que a análise: ela é descartada, não propagada.

    Além do JSON quebrado, o conteúdo pode ser válido e não ser um objeto — uma lista ou ``null``
    passam pelo ``json.loads`` e só estouram depois, no ``state.get``, já dentro da restauração.
    """
    (ctx.root / backend.STEP_STATE_FILENAME).write_text(content, encoding="utf-8")
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")

    assert backend.read_step_state(ctx) == {}
    assert backend.restore_analysis(ctx).references[0].endswith("ref_0_a.fasta")
    # e gravar de novo recompõe o arquivo em vez de estourar
    backend.record_step_state(ctx, ref_swap=True)
    assert backend.read_step_state(ctx) == {"ref_swap": True}


def test_a_step_state_that_cannot_be_written_does_not_break_the_step(ctx, warned, monkeypatch):
    """
    Gravar a seleção falha em silêncio, como o histórico.

    Os arquivos já estão salvos e a etapa está na tela; perder a lembrança de uma escolha custa um
    reenvio na próxima sessão, derrubar a etapa custa o trabalho desta.
    """
    real_write = Path.write_text

    def refuse(self, *args, **kwargs):
        if self.name == backend.STEP_STATE_FILENAME:
            raise OSError("disco cheio")
        return real_write(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", refuse)

    backend.record_step_state(ctx, ref_swap=True)

    assert not (ctx.root / backend.STEP_STATE_FILENAME).exists()
    assert any(backend.STEP_STATE_FILENAME in line for line in warned)


# --------------------------------------------- a seleção vence as chegadas, que são outra conta
def test_only_the_file_the_step_still_holds_comes_back(ctx):
    """
    A causa raiz: o uploader de cada método SUBSTITUI a lista do método.

    Uma correção enviada por cima da primeira versão deixa duas chegadas e uma escolha só.
    Reconstruída das chegadas, a etapa 2 oferecia ao POA1 o arquivo abandonado ao lado do mantido,
    e os epítopos da versão descartada voltavam à contagem sem uma palavra.
    """
    _upload(ctx, ctx.inputs_dir, "bepipred3_0_v1.fasta", ">NS1_DENV1_ref\nmktAYIA\n",
            role="prediction", method="b3")
    kept = _upload(ctx, ctx.inputs_dir, "bepipred3_0_v2.fasta", ">NS1_DENV1_ref\nMKTayia\n",
                   role="prediction", method="b3")
    backend.record_step_state(ctx, sources={"b3": [kept]})

    restored = backend.restore_analysis(ctx)

    assert restored.prediction_files == {"b3": [kept]}


def test_a_step_cleared_by_hand_comes_back_cleared(ctx):
    """
    "Limpar espécies importadas" não faz arquivo nenhum chegar: só esvazia a etapa.

    É por isso que a seleção é gravada mesmo vazia. Com a conta das chegadas sozinha, reabrir
    ressuscitava exatamente o que a pessoa tinha acabado de tirar da etapa, e o botão de limpar
    durava só até a próxima sessão.
    """
    _upload(ctx, ctx.inputs_dir, "bepipred3_0_v1.fasta", ">NS1_DENV1_ref\nmktAYIA\n",
            role="prediction", method="b3")
    backend.record_step_state(ctx, sources={})

    assert backend.restore_analysis(ctx).prediction_files == {}


def test_the_bepipred_version_the_automated_run_discarded_does_not_come_back(ctx):
    """
    BepiPred-2.0 e 3.0 se expulsam: a corrida automática do 3.0 descarta o -b2 da sessão.

    Nada chega e nada é apagado do disco — a mudança é só de seleção —, então o JSON do 2.0
    continuava conforme e voltava à etapa 2 junto com o FASTA do 3.0, numa análise que a pessoa
    tinha decidido rodar com um só.
    """
    backend.save_upload(ctx.inputs_dir, "raw_ok.json", _bepipred2_json(NS1_DENV1_ref="MKTAYIA"),
                        ctx=ctx, role="bepipred2_raw", original_name="ok.json")
    b3 = _upload(ctx, ctx.inputs_dir, "bepipred3.fasta", ">NS1_DENV1_ref\nmktAYIA\n",
                 role="prediction", method="b3")
    backend.record_step_state(ctx, sources={"b3": [b3]})

    restored = backend.restore_analysis(ctx)

    assert restored.prediction_files == {"b3": [b3]}


def test_the_adapted_imports_the_person_cleared_do_not_come_back(ctx):
    """
    O mesmo botão de limpar, do lado dos pares adaptados.

    Os dois arquivos do adaptador continuam em ``inputs/`` depois de limpar — nada os apaga —, e o
    par é encontrado pelo radical compartilhado, com registro ou sem ele. Só a seleção sabe que a
    pessoa não quer mais aquele JSON na etapa 2.
    """
    _step_two_import(ctx, "b.json", _bepipred2_json(Sequence="MKTAYIA"), specie="DENV1")
    backend.record_step_state(ctx, prepared_b2=[])

    restored = backend.restore_analysis(ctx)

    assert restored.prepared_b2 == []
    assert restored.references == [], "o -f reconstruído continua fora da etapa 1"


def test_the_adapted_import_the_person_kept_is_the_only_one_back(ctx):
    """Limpar e reimportar deixa dois pares no disco; a etapa 2 guarda um."""
    _step_two_import(ctx, "dv1.json", _bepipred2_json(Sequence="MKTAYIA"), specie="DENV1")
    kept = _step_two_import(ctx, "chikv.json", _bepipred2_json(Sequence="QRQMKTA"),
                            specie="CHIKV", protein="E1")
    backend.record_step_state(ctx, prepared_b2=[kept.b2_json])

    restored = backend.restore_analysis(ctx)

    assert [pr.b2_json for pr in restored.prepared_b2] == [kept.b2_json]


def test_the_order_the_person_chose_for_the_references_survives(ctx):
    """
    A ordem do -f é científica, não cosmética: ``merge_fastas`` guarda o PRIMEIRO cabeçalho igual.

    Duas versões da mesma proteína sob o mesmo cabeçalho fazem a ordem decidir qual sequência o
    POA1 analisa. As chegadas voltam em ordem cronológica, que é a ordem de quem enviou — e a
    etapa 1 oferece a ordem de quem escolheu.
    """
    a = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")
    b = _upload(ctx, ctx.inputs_dir, "ref_1_b.fasta", ">NS1_DENV2_ref\nQRQMKTA\n", role="reference")
    c = _upload(ctx, ctx.inputs_dir, "ref_2_c.fasta", ">NS1_DENV3_ref\nAYIAKQR\n", role="reference")
    backend.record_step_state(ctx, references=[c, a, b])

    assert backend.restore_analysis(ctx).references == [c, a, b]


def test_a_selected_file_gone_from_disk_takes_its_method_with_it(ctx):
    """
    Quem apaga o arquivo escolhido não quer de volta o que ele substituiu.

    A seleção decide quem é candidato e o disco decide quem sobrevive, nessa ordem: a etapa reabre
    vazia, pedindo um envio novo, em vez de oferecer ao POA1 a versão que a pessoa já tinha
    trocado e que continua em ``inputs/``.
    """
    _upload(ctx, ctx.inputs_dir, "bepipred3_0_v1.fasta", ">NS1_DENV1_ref\nmktAYIA\n",
            role="prediction", method="b3")
    kept = _upload(ctx, ctx.inputs_dir, "bepipred3_0_v2.fasta", ">NS1_DENV1_ref\nMKTayia\n",
                   role="prediction", method="b3")
    backend.record_step_state(ctx, sources={"b3": [kept]})
    Path(kept).unlink()

    assert backend.restore_analysis(ctx).prediction_files == {}


def test_mhcii_is_the_folder_whatever_the_selection_says(ctx):
    """
    O MHC-II é a única entrada que o POA1 recebe como PASTA de ``Proteína_Espécie.html``.

    Cada HTML é registrado por si — isso é a procedência —, mas a etapa não escolhe entre eles:
    entrega ``mhcii/`` inteira sempre que houver algo lá. Narrar o MHC-II pela seleção de arquivos
    faria a etapa 2 devolver a pasta vazia, e o POA1 rodaria sem a predição de célula T auxiliar.
    """
    _upload(ctx, ctx.mhcii_dir, "NS1_DENV1.html", "<html></html>", role="prediction", method="m")
    _upload(ctx, ctx.inputs_dir, "bepipred3_0_v1.fasta", ">NS1_DENV1_ref\nmktAYIA\n",
            role="prediction", method="b3")
    kept = _upload(ctx, ctx.inputs_dir, "bepipred3_0_v2.fasta", ">NS1_DENV1_ref\nMKTayia\n",
                   role="prediction", method="b3")
    backend.record_step_state(ctx, sources={"b3": [kept]})

    restored = backend.restore_analysis(ctx)

    assert restored.prediction_files == {"b3": [kept], "m": [str(ctx.mhcii_dir)]}


def test_an_analysis_older_than_the_step_state_keeps_what_the_files_gave_it(ctx):
    """
    Só os campos gravados são aplicados; o resto continua valendo pelos arquivos.

    Uma análise anterior a esta camada — ou que só renderizou a etapa 5 — tem um ``step_state``
    parcial. Tratar o que falta como "seleção vazia" esvaziaria as etapas 1 e 2 de toda análise
    feita até aqui, que é o oposto de restaurar.
    """
    ref = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
                  role="reference")
    b3 = _upload(ctx, ctx.inputs_dir, "bepipred3_0_x.fasta", ">NS1_DENV1_ref\nmktAYIA\n",
                 role="prediction", method="b3")
    backend.record_step_state(ctx, params_poa2={"t": 90, "objective": "conserved"})

    restored = backend.restore_analysis(ctx)

    assert restored.references == [ref]
    assert restored.prediction_files == {"b3": [b3]}
    assert restored.params == {"params_poa2": {"t": 90, "objective": "conserved"}}


# --------------------------------------------------------- a caixa de inverter e os formulários
def test_the_swapped_headers_box_comes_back_ticked(ctx):
    """
    Com a caixa marcada, o -f da etapa 1 é ``proteins_normalized.fasta``, não os envios crus.

    Reabrir com a caixa desmarcada mostrava os mesmos arquivos e reconstruía o -f sem inverter:
    o POA1 passava a ler ``DENV1`` como proteína e ``NS1`` como espécie, e o resultado saía com
    outro conjunto de espécies sem que nada na tela tivesse mudado.
    """
    backend.record_step_state(ctx, ref_swap=True)

    assert backend.restore_analysis(ctx).ref_swap is True


def test_a_box_never_ticked_comes_back_unticked(ctx):
    """O padrão é não inverter: ausência de registro não pode virar inversão."""
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")

    assert backend.restore_analysis(ctx).ref_swap is False


def test_a_ref_swap_that_is_not_a_boolean_is_ignored(ctx):
    """
    O arquivo é editável à mão e vem de outra máquina: o tipo é conferido antes de virar estado.

    ``"sim"`` é verdadeiro em Python, então aceitar qualquer coisa inverteria cabeçalhos por causa
    de um valor que o código nunca escreveu.
    """
    (ctx.root / backend.STEP_STATE_FILENAME).write_text('{"ref_swap": "sim"}', encoding="utf-8")

    assert backend.restore_analysis(ctx).ref_swap is False


def test_the_forms_come_back_beside_the_files_they_produced(ctx):
    """
    Os parâmetros são estado da análise tanto quanto os arquivos dela.

    Os HTML do MHC-II saíram de alelos específicos e a planilha do POA2 de filtros específicos.
    Devolver os arquivos ao lado dos padrões iniciais convida a repetir a etapa com entradas que
    não correspondem ao que está na tela — e, no POA2, a sobrescrever a planilha com outra.
    """
    backend.record_step_state(ctx, params_poa1={"bmin": 8, "bmax": 20, "e": "y"},
                              params_mhcii={"alleles": "HLA-DRB1*01:01"})
    backend.record_step_state(ctx, params_poa2={"t": 90, "objective": "unique", "rf": 3})

    params = backend.restore_analysis(ctx).params

    assert params["params_poa1"] == {"bmin": 8, "bmax": 20, "e": "y"}
    assert params["params_mhcii"] == {"alleles": "HLA-DRB1*01:01"}
    assert params["params_poa2"] == {"t": 90, "objective": "unique", "rf": 3}


def test_only_the_form_keys_are_taken_as_params(ctx):
    """
    ``params`` é copiado para o estado da sessão chave a chave, então só o que é formulário entra.

    O arquivo guarda seleções de arquivo no mesmo nível ('references', 'sources'); varrê-lo
    inteiro escreveria uma lista de caminhos por cima de um campo de formulário da interface.
    """
    (ctx.root / backend.STEP_STATE_FILENAME).write_text(
        json.dumps({"references": ["inputs/a.fasta"], "ref_swap": True,
                    "params_poa1": {"bmin": 8}, "params_quebrado": "não é formulário"}),
        encoding="utf-8")

    assert backend.restore_analysis(ctx).params == {"params_poa1": {"bmin": 8}}


def test_an_analysis_with_nothing_recorded_restores_no_params(ctx):
    """Sem registro, os formulários ficam como a interface os inicializa."""
    assert backend.restore_analysis(ctx).params == {}


def test_a_remembered_form_alone_is_not_something_to_restore(ctx):
    """
    ``anything`` conta arquivos, e é ele que decide se a reabertura tem o que dizer.

    Uma pasta que só guarda preferências não devolve arquivo nenhum às etapas: anunciá-la como
    "itens recarregados" mostraria um resumo vazio, e aplicar os formulários sem os arquivos que
    os justificam deixaria a tela contando uma análise que não existe.
    """
    backend.record_step_state(ctx, ref_swap=True, params_poa1={"bmin": 8})

    restored = backend.restore_analysis(ctx)

    assert restored.ref_swap is True and restored.params
    assert not restored.anything
    assert restored.summary() == []


# --------------------------------------------------------------- o resultado do POA2 lido do disco
def _poa2_spreadsheet(ctx, criterion, *, rows=2, rf=None, mtime=None):
    """
    A planilha que ``run_poa2`` deixa em ``poa2_out/``, com o critério no nome.

    O operador não está dentro do arquivo: ``>`` e ``<`` não são legais em nome de arquivo no
    Windows, então o POA2 escreve 'ge'/'lt' — e é do nome que a reabertura tem de tirá-lo de volta.
    """
    xlsx = ctx.poa2_out_dir / f"POA2_analysis_{criterion}.xlsx"
    pd.DataFrame({"Specie": ["DENV1"] * rows,
                  "Peptide Sequence": ["MKTAYIA"] * rows}).to_excel(xlsx, index=False)
    if rf is not None:
        (ctx.poa2_out_dir / f"{xlsx.stem}_rf{rf}.fasta").write_text(
            ">NS1_DENV1_B_1_7\nMKTAYIA\n", encoding="utf-8")
    if mtime is not None:
        os.utime(xlsx, (mtime, mtime))
    return xlsx


def test_a_finished_analysis_comes_back_with_its_poa2_result(ctx):
    """
    Reabrir uma análise concluída trazia a etapa 3 de volta e deixava a 6 dizendo que o POA2 não
    tinha rodado, com a planilha que ela pede parada na pasta.

    Rodar de novo para ler um arquivo que já está lá é o trabalho que a reabertura existe para
    poupar — e, pior, uma segunda corrida pode produzir uma planilha DIFERENTE da que está ao lado.
    """
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")
    (ctx.conservancy_csv_dir / "DENV1.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    xlsx = _poa2_spreadsheet(ctx, "ge70", rows=2, rf=1)

    restored = backend.restore_analysis(ctx)

    assert restored.poa2_result.xlsx_path == str(xlsx)
    assert len(restored.poa2_result.results) == 2
    # a planilha é lida de verdade, não reconstruída do nome
    assert list(restored.poa2_result.results.columns) == ["Specie", "Peptide Sequence"]
    assert restored.poa2_result.fasta_path.endswith("POA2_analysis_ge70_rf1.fasta")
    assert "Etapa 5: 2 epítopo(s) selecionado(s) (>=)" in restored.summary()


@pytest.mark.parametrize("criterion, symbol", [("ge70", ">="), ("lt70", "<")])
def test_the_operator_in_the_file_name_comes_back_with_the_spreadsheet(ctx, criterion, symbol):
    """
    ``>=`` e ``<`` são perguntas opostas: epítopos conservados contra epítopos únicos.

    O operador não sobrevive dentro do .xlsx, só no nome — e é ele que a etapa 6 mostra. Recuperar
    a planilha sem o operador rotularia um conjunto de únicos como conservado.
    """
    _poa2_spreadsheet(ctx, criterion)

    assert backend.poa2_result_from_disk(ctx).type_symbol == symbol


def test_the_most_recent_spreadsheet_is_the_one_the_step_shows(ctx):
    """
    Rodar a etapa 5 com o outro operador deixa as duas planilhas na pasta.

    Elas têm nomes diferentes justamente para que a segunda não apague a primeira; a etapa, porém,
    mostra uma só, e tem de ser a última rodada — senão reabrir volta no tempo e mostra um
    resultado que a pessoa já substituiu.
    """
    _poa2_spreadsheet(ctx, "ge70", rows=2, rf=1, mtime=1_000_000)
    _poa2_spreadsheet(ctx, "lt70", rows=1, rf=3, mtime=2_000_000)

    result = backend.poa2_result_from_disk(ctx)

    assert result.type_symbol == "<" and len(result.results) == 1
    # e o -rf que volta é o da planilha escolhida, não o primeiro da pasta
    assert result.fasta_path.endswith("POA2_analysis_lt70_rf3.fasta")


def test_a_spreadsheet_without_its_fasta_still_comes_back(ctx):
    """O ``-rf`` é opcional na etapa 5, e sua ausência não pode cancelar o resultado inteiro."""
    _poa2_spreadsheet(ctx, "ge70")

    result = backend.poa2_result_from_disk(ctx)

    assert result is not None and result.fasta_path is None


def test_an_analysis_that_never_ran_the_poa2_restores_no_result(ctx):
    """Sem planilha, a etapa 6 tem de continuar dizendo que falta rodar o POA2."""
    assert backend.poa2_result_from_disk(ctx) is None
    assert backend.restore_analysis(ctx).poa2_result is None
    # um .xlsx qualquer na pasta não é a planilha do POA2
    (ctx.poa2_out_dir / "conferencia.xlsx").write_text("x", encoding="utf-8")
    assert backend.poa2_result_from_disk(ctx) is None


def test_an_unreadable_spreadsheet_does_not_break_the_reopening(ctx, warned):
    """
    Nenhum arquivo que a reabertura apenas esperava encontrar pode derrubá-la.

    Uma planilha truncada ou pela metade faz o ``pandas`` levantar — e não um ``OSError``, mas o
    que o motor de leitura decidir —, então a análise inteira deixaria de abrir por causa de um
    arquivo que nem sequer é entrada de etapa nenhuma.
    """
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")
    (ctx.poa2_out_dir / "POA2_analysis_ge70.xlsx").write_text("isto não é uma planilha",
                                                              encoding="utf-8")

    restored = backend.restore_analysis(ctx)

    assert restored.poa2_result is None
    assert restored.references[0].endswith("ref_0_a.fasta")
    assert any("POA2_analysis_ge70.xlsx" in line for line in warned)


def test_anything_counts_what_the_summary_reports(ctx):
    """
    ``summary()`` anunciava a etapa 5 e ``anything`` a ignorava.

    ``_restore_analysis_files`` desiste quando ``anything`` é falso, então uma pasta cujo único
    artefato sobrevivente é a planilha do POA2 produzia um resumo com "Etapa 5: ..." que ninguém
    chegava a ver, e a etapa 6 continuava dizendo que o POA2 não rodou — exatamente o defeito que
    ``poa2_result_from_disk`` existe para corrigir, de pé no caso em que as entradas não estão
    mais lá. O mesmo valia para ``poa1_result``.
    """
    _poa2_spreadsheet(ctx, "ge70")

    restored = backend.restore_analysis(ctx)

    assert any("Etapa 5" in line for line in restored.summary())
    assert restored.anything


# ----------------------------------------------- o que o POA1 monta não é entrada de ninguém
@pytest.mark.parametrize("merged", sorted({name for _, name in backend._MERGERS.values()}))
def test_every_name_consolidate_sources_writes_is_an_assembled_input(merged):
    """
    Os nomes são lidos da própria tabela que os escreve, para que as duas pontas não se separem.

    Acrescentar um método a ``_MERGERS`` sem ensinar a reabertura a reconhecer o arquivo montado
    faria o POA1 receber as partes mais a soma delas, contando cada epítopo duas vezes.
    """
    assert backend._is_assembled_input(merged)


@pytest.mark.parametrize("name", [
    "bepipred3_0_DENV1.fasta",     # uma parte, que carrega o mesmo radical da soma
    "others_0_x.fasta",
    "proteins_all.fasta",          # também é uma união, mas da etapa 1: é o -f que o POA1 recebe
    "ref_0_minhas.fasta",
])
def test_a_file_somebody_supplied_is_not_an_assembled_input(name):
    """
    Descartar de menos conta epítopos duas vezes; descartar de mais perde a etapa inteira.

    ``proteins_all.fasta`` é a armadilha do meio: ele é montado pela interface, como os arquivos
    de ``_MERGERS``, mas é entrada de verdade — o ``-f`` — e some da etapa 1 se for filtrado junto.
    """
    assert not backend._is_assembled_input(name)


def test_the_merge_is_not_offered_beside_the_parts_it_was_made_of(ctx):
    """
    'others_merged.fasta' começa com o radical 'others' e termina em '.fasta'.

    Numa análise anterior ao registro, portanto, ele se classifica como predição de outros
    preditores ao lado das próprias partes — e, se o radical não casasse, cairia no fallback de
    FASTA e entraria nas proteínas de referência. As duas leituras somam o mesmo arquivo duas
    vezes à análise.
    """
    for i, sp in enumerate(["DENV1", "DENV2"]):
        (ctx.inputs_dir / f"others_{i}_{sp}.fasta").write_text(
            f">NS1_{sp}_OTHER_1_7\nMKTAYIA\n", encoding="utf-8")
    parts = [str(p) for p in sorted(ctx.inputs_dir.glob("others_?_*.fasta"))]
    backend.consolidate_sources({"x": parts}, ctx.inputs_dir)
    assert (ctx.inputs_dir / "others_merged.fasta").is_file()

    restored = backend.restore_analysis(ctx)

    assert restored.prediction_files["x"] == parts
    assert restored.references == []


# ------------------------------------------- quando um JSON do BepiPred-2.0 é o próprio -b2
@pytest.mark.parametrize("name, antigens, usable", [
    ("conforme.json", {"NS1_DENV1_ref": "MKTAYIA"}, True),
    ("generico.json", {"Sequence": "MKTAYIA"}, False),      # o que o site devolve por padrão
    ("so_especie.json", {"denv1": "MKTAYIA"}, False),
    ("dois_campos.json", {"E_DENV1": "MKTAYIA"}, False),    # falta o ID, e o parser pede três
    ("misto.json", {"NS1_DENV1_ref": "MKTAYIA", "Sequence": "QRQMKTA"}, False),
    ("sem_antigenos.json", {}, False),
])
def test_only_keys_in_protein_specie_id_make_a_json_the_b2_input(tmp_path, name, antigens, usable):
    """
    Devolver um JSON não mapeado seria pior que não devolver nada: os epítopos caem sob ``nan``.

    Três casos que uma checagem mais frouxa deixaria passar. O misto, porque basta uma chave
    genérica para que aquela sequência perca a espécie — a conformidade é do arquivo, não da
    maioria das chaves. O de dois campos, porque ``Proteína_Espécie`` parece conforme e não é. E o
    sem antígenos, porque ele não tem chave NENHUMA fora do padrão: perguntar só "sobrou alguma
    chave ruim?" aceitaria um arquivo vazio e entregaria ao POA1 um ``-b2`` sem uma predição.
    """
    path = tmp_path / name
    path.write_bytes(_bepipred2_json(**antigens))

    assert backend._usable_as_b2(path) is usable


@pytest.mark.parametrize("name, content", [
    ("quebrado.json", "{não é json"),            # ValueError na leitura
    ("lista.json", '[{"antigens": {}}]'),        # AttributeError: list não tem .get
    ("nulo.json", '{"antigens": null}'),         # TypeError: None não é iterável
])
def test_a_malformed_bepipred2_json_only_warns(tmp_path, warned, name, content):
    """
    A resposta é uma propriedade do arquivo, então é relida a cada reabertura — inclusive em
    arquivos que ninguém mais abre.

    Nenhuma dessas três exceções é ``OSError``, ``ValueError`` ou ``KeyError`` ao mesmo tempo, e
    um ``except`` estreito faria um arquivo ilegível em ``inputs/`` impedir a análise inteira de
    abrir. O arquivo sai da etapa 2 com um aviso que o nomeia, que é o que permite achá-lo.
    """
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")

    assert backend._usable_as_b2(path) is False
    assert any(name in line and "BepiPred-2.0" in line for line in warned)
