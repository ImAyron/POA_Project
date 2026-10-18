"""A seleção atual de cada etapa, gravada em ``step_state.json``.

``input_history.json`` registra CHEGADAS: que um arquivo entrou na análise, e como o quê. É um
fato sobre o arquivo que nunca deixa de ser verdade. Só que a interface trabalha por SELEÇÃO: o
uploader da etapa 2 substitui a lista do seu método, BepiPred-2.0 e 3.0 se expulsam porque são um
método B só, e "Limpar espécies importadas" descarta as importações adaptadas sem que arquivo
nenhum chegue. Nenhuma dessas escolhas deixa rastro num registro de chegadas.

Por isso uma análise reconstruída só das chegadas ressuscitava o que a pessoa tinha abandonado — o
JSON que ela substituiu, a importação que a corrida automática do BepiPred-3.0 descartou — e
entregava tudo junto ao POA1, somando epítopos de arquivos que não estavam mais em etapa nenhuma.
``step_state.json`` é a conta que grava a escolha; é ela que estes casos travam.
"""
import json
from pathlib import Path

import pytest

from poa.gui import backend


@pytest.fixture
def ctx(tmp_path, monkeypatch):
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    return backend.create_run(tmp_path, "teste step state")


def _upload(ctx, dest, name, text, **provenance):
    return backend.save_upload(dest, name, text.encode("utf-8"), ctx=ctx, **provenance)


def _bepipred2_json(**antigens):
    """Um export do BepiPred-2.0, uma chave de antígeno por sequência."""
    return json.dumps({
        "info": {"failedjobs": 0, "size": len(antigens)},
        "antigens": {key: {"AA": list(seq), "BepiPred-2.0 score": [0.1] * len(seq)}
                     for key, seq in antigens.items()},
    }).encode("utf-8")


# --------------------------------------------------------- a seleção substitui; a chegada acumula
def test_the_current_selection_replaces_the_previous_one_of_the_same_step(ctx):
    """
    O ponto inteiro da camada: o histórico acumula, o step_state substitui.

    O uploader da etapa 1 troca ``ss.ref_raw`` em vez de somar a ele. Os dois FASTAs chegaram, e as
    duas chegadas continuam registradas — mas só o segundo está selecionado, e devolver o primeiro
    junto mesclaria em proteins_all.fasta a sequência que a pessoa reabriu a análise para trocar.
    """
    velho = _upload(ctx, ctx.inputs_dir, "ref_0_velho.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
                    role="reference")
    novo = _upload(ctx, ctx.inputs_dir, "ref_0_novo.fasta", ">NS1_DENV1_ref\nQRQMKTA\n",
                   role="reference")
    backend.record_step_state(ctx, references=[velho])
    backend.record_step_state(ctx, references=[novo])

    restored = backend.restore_analysis(ctx)

    assert [Path(p).name for p in restored.references] == ["ref_0_novo.fasta"]
    # o arquivo abandonado continua no disco e na conta das chegadas: quem mudou foi a seleção
    assert Path(velho).is_file()
    assert len(backend.read_input_history(ctx)) == 2


def test_a_method_the_step_evicted_is_not_offered_to_poa1_again(ctx):
    """
    BepiPred-2.0 e 3.0 se expulsam — ``_automate_bepipred3`` faz ``poa1_sources.pop("b2")``.

    Sem a seleção gravada, as duas chegadas voltavam como dois métodos prontos e o POA1 contava os
    epítopos de células B duas vezes, a partir de um arquivo que a etapa 2 já tinha descartado.
    """
    _upload(ctx, ctx.inputs_dir, "raw_ok.json",
            _bepipred2_json(NS1_DENV1_ref="MKTAYIA").decode("utf-8"),
            role="bepipred2_raw", original_name="ok.json")
    b3 = _upload(ctx, ctx.inputs_dir, "bepipred3.fasta", ">NS1_DENV1_ref\nmktAYIa\n",
                 role="prediction", method="b3")
    assert sorted(backend.restore_analysis(ctx).prediction_files) == ["b2", "b3"]

    backend.record_step_state(ctx, sources={"b3": [b3]})

    assert sorted(backend.restore_analysis(ctx).prediction_files) == ["b3"]


def test_clearing_the_imported_species_is_recorded_as_an_empty_selection(ctx):
    """
    "Limpar espécies importadas" descarta as adaptações sem que arquivo nenhum chegue ou saia.

    É a escolha que menos rastro deixa: os dois arquivos do adaptador continuam em ``inputs/`` com
    a procedência gravada, então só a seleção pode dizer que a etapa 2 não os tem mais.
    """
    prep = backend.import_bepipred2(ctx, "b.json", _bepipred2_json(Sequence="MKTAYIA"),
                                    "DENV1", protein="E")
    assert len(backend.restore_analysis(ctx).prepared_b2) == 1

    backend.record_step_state(ctx, prepared_b2=[])

    assert backend.restore_analysis(ctx).prepared_b2 == []
    assert Path(prep.b2_json).is_file(), "limpar a etapa não apaga o arquivo"


def test_several_species_of_one_method_come_back_side_by_side(ctx):
    """
    "Substituir" é a lista inteira do método, nunca "só o último sobrevive".

    ``poa1_sources`` guarda lista justamente porque uma corrida multi-espécie produz um arquivo por
    vírus, e eles são mesclados num só antes de ir ao POA1. Aqui três espécies chegaram e o
    uploader ficou com duas: se a seleção guardasse um caminho só, toda análise multi-espécie se
    reduziria à última espécie enviada — o erro que o merge por método existe para evitar.
    """
    partes = [_upload(ctx, ctx.inputs_dir, f"bepipred3_{i}_{sp}.fasta", f">NS1_{sp}_ref\nMKTAYIA\n",
                      role="prediction", method="b3")
              for i, sp in enumerate(["DENV1", "DENV2", "DENV3"])]

    backend.record_step_state(ctx, sources={"b3": partes[:2]})

    restored = backend.restore_analysis(ctx)
    assert [Path(p).name for p in restored.prediction_files["b3"]] == [
        "bepipred3_0_DENV1.fasta", "bepipred3_1_DENV2.fasta"]


def test_a_step_whose_selection_is_empty_restores_nothing_for_it(ctx):
    """
    Uma etapa sem seleção nenhuma é uma resposta, não a ausência de uma.

    ``step_predictions`` grava a seleção mesmo vazia, que é exatamente o que limpar a etapa produz;
    tratar o vazio como "nada foi gravado" devolveria os arquivos que a pessoa acabou de tirar.
    """
    _upload(ctx, ctx.inputs_dir, "others_0_x.fasta", ">NS1_DENV1_O_1_3\nMKT\n",
            role="prediction", method="x")

    backend.record_step_state(ctx, sources={})

    assert backend.restore_analysis(ctx).prediction_files == {}


def test_recording_one_step_does_not_narrow_the_steps_that_said_nothing(ctx):
    """
    Cada etapa grava o seu campo; o que ninguém gravou continua vindo dos arquivos.

    É o que mantém as análises anteriores a esta camada abrindo como antes — e, dentro de uma
    sessão, o que impede a etapa 2 de apagar a etapa 1 ao gravar a própria seleção.
    """
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")
    _upload(ctx, ctx.inputs_dir, "ref_0_b.fasta", ">NS2_DENV1_ref\nQRQMKTA\n", role="reference")
    x = _upload(ctx, ctx.inputs_dir, "others_0_x.fasta", ">NS1_DENV1_O_1_3\nMKT\n",
                role="prediction", method="x")
    _upload(ctx, ctx.inputs_dir, "others_1_y.fasta", ">NS1_DENV1_O_4_6\nAYI\n",
            role="prediction", method="x")

    backend.record_step_state(ctx, sources={"x": [x]})

    restored = backend.restore_analysis(ctx)
    assert "references" not in backend.read_step_state(ctx)
    # a etapa 2 gravou e foi estreitada; a etapa 1 não gravou e continua com as duas chegadas
    assert [Path(p).name for p in restored.prediction_files["x"]] == ["others_0_x.fasta"]
    assert [Path(p).name for p in restored.references] == ["ref_0_a.fasta", "ref_0_b.fasta"]


def test_a_selected_file_that_is_gone_from_disk_is_dropped(ctx):
    """
    A seleção diz o que a etapa tinha; o disco decide o que ela ainda tem.

    Devolver um caminho morto seria pior que devolver menos: a etapa 1 lê cada FASTA para montar o
    ``-f``, e o arquivo apagado entre duas sessões derrubaria a reabertura inteira.
    """
    a, b, c = [_upload(ctx, ctx.inputs_dir, f"ref_0_{n}.fasta", f">NS1_DENV{i}_ref\nMKTAYIA\n",
                       role="reference") for i, n in enumerate("abc")]
    backend.record_step_state(ctx, references=[a, b])
    Path(a).unlink()

    restored = backend.restore_analysis(ctx)

    assert [Path(p).name for p in restored.references] == ["ref_0_b.fasta"]
    assert Path(c).is_file(), "o terceiro continua no disco; o que o tirou da etapa foi a seleção"


def test_mhcii_is_the_folder_itself_and_survives_a_selection_that_omits_it(ctx):
    """
    MHC-II é a pasta ``mhcii/``, não uma escolha entre arquivos.

    A etapa entrega a pasta ao POA1 sempre que ela tem algo dentro, então o ``m`` da seleção é o
    caminho da própria pasta e não passa pelo teste de "é um arquivo" — sem a exceção, gravar a
    seleção da etapa 2 derrubaria o método que ela mesma acabou de preencher.
    """
    _upload(ctx, ctx.mhcii_dir, "NS1_DENV1.html", "<html></html>",
            role="prediction", method="m")
    b3 = _upload(ctx, ctx.inputs_dir, "bepipred3_0_a.fasta", ">NS1_DENV1_ref\nmktAYIa\n",
                 role="prediction", method="b3")
    _upload(ctx, ctx.inputs_dir, "bepipred3_1_b.fasta", ">NS2_DENV1_ref\nqrqMKTa\n",
            role="prediction", method="b3")

    backend.record_step_state(ctx, sources={"b3": [b3], "m": [str(ctx.mhcii_dir)]})

    restored = backend.restore_analysis(ctx)
    assert restored.prediction_files["m"] == [str(ctx.mhcii_dir)]
    assert [Path(p).name for p in restored.prediction_files["b3"]] == ["bepipred3_0_a.fasta"]


# --------------------------------------------------------------------------------- os caminhos
def test_paths_are_recorded_relative_to_the_analysis_folder(ctx):
    """
    A mesma pasta é aberta do Windows e da WSL, onde o caminho absoluto é outro.

    O registro das chegadas já grava relativo pelo mesmo motivo; a seleção grava os mesmos
    arquivos, e um absoluto aqui tornaria a etapa 2 abrir vazia do outro lado do vai-e-vem entre os
    dois ambientes, que é o fluxo de trabalho do projeto.
    """
    path = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
                   role="reference")

    backend.record_step_state(ctx, references=[path])

    gravado = json.loads((ctx.root / backend.STEP_STATE_FILENAME).read_text(encoding="utf-8"))
    assert gravado["references"] == ["inputs/ref_0_a.fasta"]
    assert ":" not in gravado["references"][0] and not gravado["references"][0].startswith("/")
    # e o caminho absoluto é remontado na leitura, que é o que a etapa recebe de volta
    assert backend.restore_analysis(ctx).references == [path]


def test_every_path_inside_a_recorded_value_is_rewritten(ctx):
    """
    A seleção não é uma lista de caminhos: é um dicionário por método, e os formulários trazem
    caminho dentro de si (o ``-rf`` do POA2).

    Relativizar só o topo gravaria absoluto tudo o que estivesse um nível abaixo, que é onde está a
    maior parte dos caminhos desta camada.
    """
    dentro = str(ctx.inputs_dir / "a.fasta")
    fora = str(Path(ctx.root).parent / "fora.fasta")

    recortado = backend._relativise(ctx, {"b3": (dentro,), "params": {"rf": dentro, "t": 70}})

    assert recortado == {"b3": ["inputs/a.fasta"], "params": {"rf": "inputs/a.fasta", "t": 70}}
    # um arquivo que não está na pasta da análise não tem forma relativa: fica como veio
    assert backend._relativise(ctx, fora) == fora


# ----------------------------------------------------------------------------------- a gravação
def test_recording_the_same_state_again_does_not_touch_the_file(ctx):
    """
    O Streamlit re-executa a etapa inteira a cada clique, e ela grava a seleção toda vez.

    Reescrever o arquivo a cada interação de uma análise onde nada mudou é escrita no disco do
    usuário por conta nenhuma — e apaga a única pista de quando a seleção mudou de fato. O conteúdo
    sozinho não distingue "não gravou" de "gravou igual"; o mtime distingue.
    """
    path = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
                   role="reference")
    backend.record_step_state(ctx, references=[path])
    target = ctx.root / backend.STEP_STATE_FILENAME
    antes = (target.stat().st_mtime_ns, target.read_bytes())

    backend.record_step_state(ctx, references=[path])

    assert (target.stat().st_mtime_ns, target.read_bytes()) == antes
    # a troca de verdade ainda passa — o que foi dispensado foi a reescrita idêntica
    backend.record_step_state(ctx, references=[])
    assert target.read_bytes() != antes[1]


def test_recording_a_step_keeps_what_the_other_steps_had_recorded(ctx):
    """
    Um arquivo só, escrito por três etapas em momentos diferentes.

    Cada chamada passa os campos da sua etapa; gravar substituindo o arquivo inteiro faria a etapa
    5 apagar a seleção da etapa 1 no instante em que alguém mexesse num filtro do POA2.
    """
    backend.record_step_state(ctx, references=["inputs/a.fasta"], ref_swap=True)
    backend.record_step_state(ctx, sources={"x": ["inputs/x.fasta"]})
    backend.record_step_state(ctx, params_poa2={"t": 70.0})

    assert backend.read_step_state(ctx) == {
        "references": ["inputs/a.fasta"],
        "ref_swap": True,
        "sources": {"x": ["inputs/x.fasta"]},
        "params_poa2": {"t": 70.0},
    }


def test_an_unreadable_step_state_is_ignored_and_the_next_write_rebuilds_it(ctx):
    """
    Um JSON truncado não pode impedir a análise de abrir.

    A seleção é uma conveniência: sem ela a reabertura ainda tem as chegadas e os arquivos em
    disco, que é o que esta camada veio estreitar — e estreitar a partir de lixo seria pior que não
    estreitar.
    """
    _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n", role="reference")
    (ctx.root / backend.STEP_STATE_FILENAME).write_text("{não é json", encoding="utf-8")

    assert backend.read_step_state(ctx) == {}
    assert [Path(p).name for p in backend.restore_analysis(ctx).references] == ["ref_0_a.fasta"]

    backend.record_step_state(ctx, ref_swap=True)

    assert backend.read_step_state(ctx) == {"ref_swap": True}


def test_a_step_state_that_is_not_an_object_is_ignored(ctx):
    """
    JSON válido não é o bastante: a camada lê campos por nome.

    Uma lista no topo passa pelo ``json.loads`` e só quebraria mais adiante, em ``state.get`` —
    dentro da reabertura, longe do arquivo que a causou.
    """
    (ctx.root / backend.STEP_STATE_FILENAME).write_text("[1, 2]", encoding="utf-8")

    assert backend.read_step_state(ctx) == {}
    assert backend.restore_analysis(ctx).prediction_files == {}


def test_a_step_state_that_cannot_be_written_never_blocks_the_upload(ctx):
    """
    Gravar a seleção é um acréscimo, nunca uma condição para salvar o arquivo.

    Mesmo motivo do histórico de chegadas: o arquivo da pessoa vale mais que a preferência que o
    acompanha, e a etapa tem de continuar funcionando com a pasta sem permissão de escrita.
    """
    (ctx.root / backend.STEP_STATE_FILENAME).mkdir()   # o caminho existe e não é gravável

    backend.record_step_state(ctx, ref_swap=True)      # não levanta

    path = _upload(ctx, ctx.inputs_dir, "ref_0_a.fasta", ">NS1_DENV1_ref\nMKTAYIA\n",
                   role="reference")
    assert Path(path).is_file()
    assert backend.read_step_state(ctx) == {}
    assert [Path(p).name for p in backend.restore_analysis(ctx).references] == ["ref_0_a.fasta"]


# ------------------------------------------------------------------- o que a etapa tinha na tela
def test_the_swap_box_comes_back_the_way_it_was_left(ctx):
    """
    A caixa de inverter os dois primeiros campos do cabeçalho decide qual arquivo é o ``-f``.

    Com ela marcada, a etapa 1 roda com proteins_normalized.fasta; reabrir com a caixa desmarcada
    ao lado desse mesmo ``-f`` mostra uma tela que não corresponde ao que a análise rodou, e
    remontar a referência dali entrega ao POA1 outro conjunto de espécies, sem dizer nada.
    """
    backend.record_step_state(ctx, ref_swap=True)
    assert backend.restore_analysis(ctx).ref_swap is True

    backend.record_step_state(ctx, ref_swap=False)
    assert backend.restore_analysis(ctx).ref_swap is False


def test_only_the_forms_come_back_as_params(ctx):
    """
    ``restored.params`` é copiado direto para as chaves ``params_*`` da sessão.

    Os HTMLs da etapa 2 foram produzidos com certos alelos e a planilha da etapa 5 com certos
    filtros: oferecer os arquivos de volta ao lado dos valores iniciais convida a re-rodar a etapa
    com entradas que não são as da tela. O que não é formulário tem campo próprio e não pode
    escorregar para cá — ``sources`` viraria ``ss.sources``, uma chave que ninguém lê.
    """
    backend.record_step_state(ctx, params_poa1={"bmin": 8, "e": "y"}, params_poa2={"t": 70.0},
                              ref_swap=True, sources={"x": []}, params_quebrado="nem dicionário")

    params = backend.restore_analysis(ctx).params

    assert params == {"params_poa1": {"bmin": 8, "e": "y"}, "params_poa2": {"t": 70.0}}


def test_a_relative_results_root_still_records_paths_relative_to_the_analysis(tmp_path, monkeypatch):
    """
    `POA_RESULTS_DIR` é override documentado do usuário e pode ser relativo.

    Quando é, os caminhos que chegam ao registro são relativos ao diretório de trabalho, não à
    pasta da análise. `_relativise` pulava esses por não serem absolutos, e o step_state ficava
    guardando um caminho relativo ao CWD; `_apply_step_state` depois resolvia como `ctx.root / p`,
    duplicando o prefixo, não achava o arquivo e descartava a seleção inteira — a etapa reabria
    vazia com o arquivo intacto no disco.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("POA_RESULTS_DIR", "resultados")
    ctx = backend.create_run(backend.default_results_root(), "raiz relativa")
    # O `-f` unido é a porta por onde um caminho relativo chega: `_refresh_reference` o monta com
    # `str(ctx.inputs_dir / "proteins_all.fasta")`, direto da raiz não resolvida. Os uploads passam
    # por `safe_child`, que resolve, e por isso nunca exercitam este caminho.
    merged = str(ctx.inputs_dir / "proteins_all.fasta")
    assert not Path(merged).is_absolute(), "o cenário exige o caminho relativo que o app produz"
    ctx.inputs_dir.mkdir(parents=True, exist_ok=True)
    Path(merged).write_text(">NS1_DENV1_ref\nMKTAYIA\n", encoding="utf-8")

    backend.record_step_state(ctx, references=[merged])

    recorded = json.loads((ctx.root / backend.STEP_STATE_FILENAME).read_text(encoding="utf-8"))
    assert recorded["references"] == ["inputs/proteins_all.fasta"]
