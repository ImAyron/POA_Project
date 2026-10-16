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
