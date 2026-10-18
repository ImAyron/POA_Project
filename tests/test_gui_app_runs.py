"""The analysis picker as the user meets it: create a run, reopen one, never mix two."""
import json

import pytest

pytest.importorskip("streamlit.testing.v1", reason="streamlit not installed")

from datetime import date  # noqa: E402
from pathlib import Path  # noqa: E402
from types import SimpleNamespace  # noqa: E402

from streamlit.testing.v1 import AppTest  # noqa: E402

from poa.core import pipeline  # noqa: E402
from poa.gui import backend  # noqa: E402
from poa.services import conservancy_client  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")


@pytest.fixture
def results_root(tmp_path, monkeypatch):
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    return tmp_path


def _start(**state):
    at = AppTest.from_file(APP, default_timeout=90)
    for key, value in state.items():
        at.session_state[key] = value
    at.run()
    return at


def test_no_analysis_selected_shows_the_picker_and_stops(results_root):
    at = _start()

    assert not at.exception
    assert any("Escolha a análise" in str(h.value) for h in at.header)
    assert not at.sidebar.radio            # the step selector is not reachable yet


def test_creating_an_analysis_binds_the_session_to_its_folder(results_root):
    at = _start()
    at.text_input[0].set_value("Teste DENV").run()
    next(b for b in at.button if "Criar" in b.label).click().run()

    assert not at.exception
    expected = results_root / backend.run_dir_name("Teste DENV")
    assert at.session_state["run_path"] == str(expected)
    assert expected.is_dir()
    assert at.sidebar.radio               # the steps are now reachable
    assert any(expected.name in str(c.value) for c in at.sidebar.code)


def test_existing_analyses_are_offered(results_root):
    backend.create_run(results_root, "Teste A", date(2026, 1, 5))
    backend.create_run(results_root, "Teste B", date(2026, 9, 27))

    at = _start()
    options = at.selectbox[0].options if at.selectbox else []
    assert len(options) == 2
    # newest first, labelled with the name the user typed and how far it got
    assert options[0].startswith("Teste B — 27/09/2026")
    assert options[1].startswith("Teste A — 05/01/2026")


def test_every_step_renders_inside_an_analysis(results_root):
    """ctx is now bound by the picker rather than at import; no step may break on that."""
    ctx = backend.create_run(results_root, "Teste A", date(2026, 9, 27))
    at = _start(run_path=str(ctx.root))
    steps = at.sidebar.radio[0].options

    assert len(steps) == 7
    for step in steps:
        at.sidebar.radio[0].set_value(step).run()
        assert not at.exception, f"{step}: {at.exception}"


def test_switching_analysis_clears_the_previous_one(results_root):
    """State points at files inside a folder; carrying it over would write into the wrong run."""
    first = backend.create_run(results_root, "Teste A", date(2026, 9, 27))
    second = backend.create_run(results_root, "Teste B", date(2026, 9, 27))

    at = _start(run_path=str(first.root), proteins_path=str(first.inputs_dir / "p.fasta"),
                conservancy_ready=True)
    assert at.session_state["proteins_path"]

    next(b for b in at.sidebar.button if "Trocar" in b.label).click().run()
    at.selectbox[0].set_value(1).run()           # index 1: the other analysis
    next(b for b in at.button if "Abrir" in b.label).click().run()

    assert at.session_state["run_path"] in (str(first.root), str(second.root))
    assert at.session_state["proteins_path"] is None
    assert at.session_state["conservancy_ready"] is False


def test_reopening_an_analysis_loads_its_files_onto_their_steps(results_root):
    """
    O pedido, de ponta a ponta: reabrir e encontrar tudo no lugar.

    Reabrir restaurava só os CSVs de conservância e a referência, e atrás de um botão na etapa 5.
    Os arquivos de predição enviados à mão na etapa 2, o conjunto de comparação da etapa 4 e a
    estrutura da etapa 7 ficavam no disco enquanto a interface mostrava campos vazios.
    """
    ctx = backend.create_run(results_root, "Teste reabrir", date(2026, 10, 16))
    backend.save_upload(ctx.inputs_dir, "ref_0_prot.fasta", b">NS1_DENV1_ref\nMKTAYIAKQRQ\n",
                        ctx=ctx, role="reference", original_name="prot.fasta")
    backend.save_upload(ctx.inputs_dir, "others_0_ep.fasta", b">NS1_DENV1_OTHER_1_5\nMKTAY\n",
                        ctx=ctx, role="prediction", method="x", original_name="ep.fasta")
    backend.save_upload(ctx.inputs_dir, "comparison_DENV1.fasta", b">NS1_DENV1_w\nMKTAYIAKQRQ\n",
                        ctx=ctx, role="comparison_set", specie="DENV1", original_name="w.fasta")
    backend.save_upload(ctx.inputs_dir, "structure_1ABC.pdb", b"ATOM\nEND\n",
                        ctx=ctx, role="structure", original_name="1abc.pdb")

    at = _start()
    next(b for b in at.button if "Abrir" in b.label).click().run()

    assert not at.exception
    assert at.session_state["run_path"] == str(ctx.root)
    # cada arquivo na etapa de onde veio
    assert at.session_state["ref_raw"][0].endswith("ref_0_prot.fasta")
    assert at.session_state["proteins_path"]
    assert at.session_state["poa1_sources"]["x"][0].endswith("others_0_ep.fasta")
    assert at.session_state["comparison_sets"]["DENV1"].endswith("comparison_DENV1.fasta")
    assert at.session_state["pdb_path"].endswith("structure_1ABC.pdb")
    # e a barra lateral diz o que voltou
    assert at.session_state["restored_source"] == "history"
    assert any("Etapa 2" in line for line in at.session_state["restored_summary"])


def test_reopening_an_old_analysis_after_a_new_upload_keeps_everything(results_root):
    """
    Uma análise anterior ao registro ganha histórico ao receber um upload qualquer.

    Daí em diante, ler só o registro escondia tudo o que já estava no disco: a reabertura mostrava
    campos vazios e o aviso de "classificado pelo nome" sumia justamente quando a restauração
    passou a fazer menos. Atinge exatamente as análises que motivaram a classificação por nome.
    """
    ctx = backend.create_run(results_root, "Teste antiga", date(2026, 10, 16))
    (ctx.inputs_dir / "ref_0_prot.fasta").write_text(">NS1_DENV1_ref\nMKTAYIAKQRQ\n",
                                                     encoding="utf-8")
    (ctx.inputs_dir / "others_0_ep.fasta").write_text(">NS1_DENV1_OTHER_1_5\nMKTAY\n",
                                                      encoding="utf-8")
    (ctx.inputs_dir / "comparison_DENV1.fasta").write_text(">NS1_DENV1_w\nMKTAYIAKQRQ\n",
                                                           encoding="utf-8")
    (ctx.inputs_dir / "structure_1ABC.pdb").write_text("ATOM\nEND\n", encoding="utf-8")
    backend.save_upload(ctx.inputs_dir, "comparison_DENV2.fasta", b">NS1_DENV2_w\nMKTAYIAKQRQ\n",
                        ctx=ctx, role="comparison_set", specie="DENV2", original_name="w2.fasta")

    at = _start()
    next(b for b in at.button if "Abrir" in b.label).click().run()

    assert not at.exception
    assert at.session_state["ref_raw"][0].endswith("ref_0_prot.fasta")
    assert at.session_state["poa1_sources"]["x"][0].endswith("others_0_ep.fasta")
    assert at.session_state["pdb_path"].endswith("structure_1ABC.pdb")
    assert sorted(at.session_state["comparison_sets"]) == ["DENV1", "DENV2"]
    # e a barra lateral diz que parte disso é palpite
    assert at.session_state["restored_source"] == "mixed"
    assert any("classificada pelo nome" in str(w.value) for w in at.sidebar.warning)


def test_the_merged_reference_is_recorded_so_reopening_stays_a_full_record(results_root):
    """
    O -f unido nasce em inputs/ sem passar por uploader nenhum, então é gravado onde nasce.

    Sem isso, a varredura o encontraria sem procedência, cairia no nome dele e pediria conferência
    numa análise cujo registro está completo — um aviso falso criado pela própria correção.
    """
    ctx = backend.create_run(results_root, "Teste duas refs", date(2026, 10, 16))
    for i, sp in enumerate(["DENV1", "DENV2"]):
        backend.save_upload(ctx.inputs_dir, f"ref_{i}_{sp}.fasta",
                            f">NS1_{sp}_ref\nMKTAYIAKQRQ\n".encode("utf-8"),
                            ctx=ctx, role="reference", original_name=f"{sp}.fasta")

    at = _start()
    next(b for b in at.button if "Abrir" in b.label).click().run()

    assert not at.exception
    assert at.session_state["proteins_path"].endswith("proteins_all.fasta")
    roles = {item["path"]: item["role"] for item in backend.read_input_history(ctx)}
    assert roles["inputs/proteins_all.fasta"] == "merged_reference"
    assert backend.restore_analysis(ctx).source == "history"
    assert not at.sidebar.warning


def test_reopening_puts_a_direct_bepipred2_json_back_on_step_2(results_root):
    """
    Um JSON do BepiPred-2.0 cujas chaves já conformam é o próprio -b2 da etapa 2.

    Gravado como "bruto" e sem ramo na restauração, ele era descartado em silêncio: uma análise
    feita pelo fluxo recomendado — submeter ao BepiPred-2.0 o FASTA que o POA0 escreve — reabria
    como se estivesse vazia.
    """
    ctx = backend.create_run(results_root, "Teste b2 direto", date(2026, 10, 16))
    raw = json.dumps({"antigens": {"NS1_DENV1_ref": {"AA": list("MKTAYIAKQRQ"),
                                                     "BepiPred-2.0 score": [0.1] * 11}}})
    backend.save_upload(ctx.inputs_dir, "raw_bepipred.json", raw.encode("utf-8"),
                        ctx=ctx, role="bepipred2_raw", original_name="bepipred.json")

    at = _start()
    next(b for b in at.button if "Abrir" in b.label).click().run()

    assert not at.exception
    assert at.session_state["poa1_sources"]["b2"][0].endswith("raw_bepipred.json")


def test_restoring_runs_once_not_on_every_interaction(results_root):
    """
    Depois da adoção, o estado da sessão é a autoridade.

    Repetir a restauração a cada reexecução desfaria o que a pessoa mexeu desde então — por
    exemplo, readicionaria um arquivo de predição que ela acabou de limpar.
    """
    ctx = backend.create_run(results_root, "Teste uma vez", date(2026, 10, 16))
    backend.save_upload(ctx.inputs_dir, "others_0_ep.fasta", b">NS1_DENV1_OTHER_1_5\nMKTAY\n",
                        ctx=ctx, role="prediction", method="x", original_name="ep.fasta")

    at = _start()
    at = next(b for b in at.button if "Abrir" in b.label).click().run()
    assert at.session_state["poa1_sources"]["x"]
    assert at.session_state["restore_pending"] is False

    # o usuário limpa a fonte e interage de novo; a restauração não pode trazê-la de volta
    at.session_state["poa1_sources"] = {}
    at = at.sidebar.radio[0].set_value("2 · Predições (auto/manual)").run()

    assert not at.exception
    assert at.session_state["poa1_sources"] == {}


# =================================================== a seleção atual de cada etapa, e seus números
# Os testes acima são sobre os arquivos que a análise contém. Estes são sobre o que cada etapa
# *escolheu* entre eles e sobre os números do formulário que acompanham essa escolha — duas coisas
# que o registro de chegadas não tem como saber: ele acumula enquanto a interface substitui
# (``ss.ref_raw`` e ``ss.poa1_sources[m]`` são trocados inteiros a cada envio), e número nenhum de
# formulário chega em ``inputs/``.

STEP2 = "2 · Predições (auto/manual)"
STEP5 = "5 · POA2 (conservação + topologia)"
STEP6 = "6 · Resultados"


def _widget(elements, label):
    return next(w for w in elements if label in w.label)


def _reopen():
    """Reabrir pela interface numa sessão nova — o app reiniciado, que é o caso real."""
    at = _start()
    assert not at.exception
    return next(b for b in at.button if "Abrir" in b.label).click().run()


def _reference(ctx, name="ref_0_prot.fasta", header="NS1_DENV1_ref"):
    return backend.save_upload(ctx.inputs_dir, name, f">{header}\nMKTAYIAKQRQ\n".encode("utf-8"),
                               ctx=ctx, role="reference", original_name=name)


def _step_one_swapped(ctx):
    """
    A etapa 1 como a pessoa a deixou: FASTA em ``Espécie_Proteína`` e a caixa de inverter marcada.

    Pela própria etapa, e não escrevendo o estado à mão, porque é ela que grava o
    ``proteins_normalized.fasta`` e a escolha que o reabrir precisa encontrar.
    """
    raw = _reference(ctx, "ref_0_denv.fasta", "DENV1_NS1_ref")
    at = _start(run_path=str(ctx.root), ref_raw=[raw])
    at = _widget(at.checkbox, "cabeçalhos").set_value(True).run()
    assert not at.exception
    assert at.session_state["proteins_path"].endswith("proteins_normalized.fasta")
    return at


def _with_conservancy(ctx):
    """A etapa 4 já feita: um epítopo, seu CSV de conservância e o limiar gravado ao lado."""
    ctx.conservancy_epitopes_dir.mkdir(parents=True, exist_ok=True)
    (ctx.conservancy_epitopes_dir / "DENV1_epitopes.fasta").write_text(
        ">DENV1_NS1_BEPIPRED_1_5\nAYIAM\n", encoding="utf-8")
    conservancy_client.run_conservancy_for_dir(
        str(ctx.conservancy_epitopes_dir),
        [("NS1_DENV1_p1", "GGAYIAMGG"), ("NS1_DENV1_p2", "TTAYIAMTT"),
         ("NS1_DENV1_p3", "GGAYIALGG")],
        threshold=70, out_dir=str(ctx.conservancy_csv_dir), operator=">=")


def _step_two_with_alleles(ctx, alleles, length):
    """A etapa 2 com os HTML do MHC-II na pasta e os alelos que os produziram no formulário."""
    ref = _reference(ctx)
    backend.save_upload(ctx.mhcii_dir, "NS1_DENV1.html", b"<html><body></body></html>",
                        ctx=ctx, role="prediction", method="m", original_name="NS1_DENV1.html")
    at = _start(run_path=str(ctx.root), ref_raw=[ref], poa1_sources={"m": [str(ctx.mhcii_dir)]})
    at = at.sidebar.radio[0].set_value(STEP2).run()
    at = _widget(at.text_input, "Alelos").set_value(alleles).run()
    at = _widget(at.number_input, "Comprimento").set_value(length).run()
    assert not at.exception
    assert at.session_state["params_mhcii"] == {"alleles": alleles, "length": length}
    return at


def test_reopening_a_swapped_analysis_gives_back_the_f_it_ran_with(results_root):
    """
    Com a inversão de cabeçalhos marcada, o `-f` da etapa 1 é o `proteins_normalized.fasta`.

    Este era o defeito bloqueante. A caixa voltava desmarcada e a etapa 1 remontava o `-f` a
    partir dos arquivos crus, cujos cabeçalhos ainda estão em `Espécie_Proteína`: a análise
    reaberta passava a ter `NS1` como espécie no lugar de `DENV1`, sem uma linha na tela, e é por
    esse campo que as etapas 3 e 4 separam as espécies.
    """
    ctx = backend.create_run(results_root, "Cabecalhos invertidos", date(2026, 10, 16))
    _step_one_swapped(ctx)

    at = _reopen()

    assert not at.exception
    assert at.session_state["ref_swap"] is True
    assert _widget(at.checkbox, "cabeçalhos").value is True
    assert at.session_state["proteins_path"].endswith("proteins_normalized.fasta")
    assert backend.species_in_fasta(at.session_state["proteins_path"]) == ["DENV1"]


def test_switching_analysis_does_not_inherit_the_header_swap(results_root):
    """
    `ref_swap` não está entre as chaves que a troca de análise limpa, e agora age sozinho.

    Enquanto reabrir deixava a etapa 1 vazia isso era inócuo: não havia o que inverter. Com os
    arquivos de volta, a caixa herdada da análise anterior reescreve os cabeçalhos da nova sem que
    ninguém tenha pedido — `NS1_DENV1_ref` vira `DENV1_NS1_ref` e a espécie lida passa a ser NS1.
    """
    swapped = backend.create_run(results_root, "Invertida", date(2026, 10, 16))
    plain = backend.create_run(results_root, "Normal", date(2026, 10, 17))   # a mais recente
    _reference(plain)

    at = _step_one_swapped(swapped)
    at = next(b for b in at.sidebar.button if "Trocar" in b.label).click().run()
    at = next(b for b in at.button if "Abrir" in b.label).click().run()

    assert not at.exception
    assert at.session_state["run_path"] == str(plain.root)
    assert at.session_state["ref_swap"] is False
    assert backend.species_in_fasta(at.session_state["proteins_path"]) == ["DENV1"]


def test_the_mhcii_parameters_come_back_with_the_files_they_produced(results_root):
    """
    Os alelos e o comprimento do peptídeo voltam junto com os HTML que saíram deles.

    Devolver os arquivos ao lado do `HLA-DRB1*01:01` inicial convida a reexecutar a etapa com
    entradas que não são as que estão na tela: os HTML na pasta foram produzidos com outros
    alelos, e nada na página diria isso.
    """
    ctx = backend.create_run(results_root, "Alelos proprios", date(2026, 10, 16))
    _step_two_with_alleles(ctx, "HLA-DRB1*04:01,HLA-DRB1*15:01", 18)

    at = _reopen()
    at = at.sidebar.radio[0].set_value(STEP2).run()

    assert not at.exception
    assert at.session_state["params_mhcii"] == {"alleles": "HLA-DRB1*04:01,HLA-DRB1*15:01",
                                                "length": 18}
    assert _widget(at.text_input, "Alelos").value == "HLA-DRB1*04:01,HLA-DRB1*15:01"
    assert _widget(at.number_input, "Comprimento").value == 18
    # o MHC-II é a pasta inteira, e é a exceção à regra da seleção: ela volta com os parâmetros
    assert at.session_state["poa1_sources"]["m"] == [str(ctx.mhcii_dir)]


def test_the_poa2_filters_come_back_with_the_csvs_they_filter(results_root):
    """
    Os filtros da etapa 5 são os parâmetros científicos da análise, e voltam com ela.

    A planilha em `poa2_out/` só é interpretável ao lado dos valores que a produziram. Reabrir
    numa `identidade mín. 60` que não foi a usada faz a etapa oferecer um reexecutar que gera uma
    planilha diferente da que está na pasta, com o mesmo nome.
    """
    ctx = backend.create_run(results_root, "Filtros proprios", date(2026, 10, 16))
    ref = _reference(ctx)
    _with_conservancy(ctx)

    at = _start(run_path=str(ctx.root), ref_raw=[ref], conservancy_ready=True)
    at = at.sidebar.radio[0].set_value(STEP5).run()
    at = _widget(at.number_input, "Identidade mín").set_value(80).run()
    at = _widget(at.number_input, "match (-m)").set_value(95).run()
    assert at.session_state["params_poa2"]["imin"] == 80

    at = _reopen()
    at = at.sidebar.radio[0].set_value(STEP5).run()

    assert not at.exception
    assert at.session_state["params_poa2"]["imin"] == 80
    assert _widget(at.number_input, "Identidade mín").value == 80
    assert _widget(at.number_input, "match (-m)").value == 95


def test_reopening_a_finished_analysis_shows_its_poa2_results_again(results_root, monkeypatch):
    """
    A planilha do POA2 está na pasta, e a etapa 6 pedia para executar de novo.

    `poa1_result_from_disk` existia sozinho: reabrir trazia a etapa 3 de volta e deixava a 6
    dizendo que o POA2 não foi executado, com a planilha que ela pede ao lado. Reexecutar para ler
    um arquivo que já está lá é o trabalho que reabrir existe para poupar — e, com os filtros da
    etapa 5 só agora restaurados, a reexecução podia gravar por cima uma planilha diferente.
    """
    ctx = backend.create_run(results_root, "POA2 concluido", date(2026, 10, 16))
    ref = _reference(ctx)
    _with_conservancy(ctx)
    # pyTMHMM é dependência nativa opcional e não decide nada do que se lê de volta.
    monkeypatch.setattr(pipeline.topology, "tmhmmAnalysis",
                        lambda args, df: df.assign(Portion_Outside=1, Portion_TM=0,
                                                   Portion_Inside=0))
    ran = pipeline.run_poa2(SimpleNamespace(
        d=str(ctx.conservancy_csv_dir), f=ref, t=70, g="y", l="", r=str(ctx.poa2_out_dir),
        imin=0, imax=100, m=0, rf=None, idf=False, strict=False))

    at = _reopen()
    at = at.sidebar.radio[0].set_value(STEP6).run()

    assert not at.exception
    assert at.session_state["poa2_result"].xlsx_path == ran.xlsx_path
    assert at.session_state["poa2_result"].type_symbol == ">="
    assert not any("Execute o POA2" in str(w.value) for w in at.main.warning)
    assert any("Etapa 5" in line for line in at.session_state["restored_summary"])


def test_reopening_gives_back_only_the_reference_that_replaced_the_other(results_root):
    """
    O uploader da etapa 1 troca `ss.ref_raw` inteiro: reenviar é dizer "não é aquele, é este".

    O registro guarda os dois, e o `-f` remontado ao reabrir ganhava uma espécie que a análise
    original não tinha. A própria etapa 5 avisa que um `-f` com espécies a mais faz um mesmo
    epítopo casar com vários registros — pior que o campo vazio que a restauração veio resolver,
    porque o resultado sai errado em vez de faltar.
    """
    ctx = backend.create_run(results_root, "Referencia trocada", date(2026, 10, 16))
    first = _reference(ctx, "ref_0_a.fasta", "NS1_DENV1_ref")
    second = _reference(ctx, "ref_0_b.fasta", "NS1_DENV2_ref")

    at = _start(run_path=str(ctx.root), ref_raw=[first])
    at.session_state["ref_raw"] = [second]       # o envio seguinte substitui a lista inteira
    at.run()
    assert not at.exception

    at = _reopen()

    assert not at.exception
    assert [Path(p).name for p in at.session_state["ref_raw"]] == ["ref_0_b.fasta"]
    assert backend.species_in_fasta(at.session_state["proteins_path"]) == ["DENV2"]


def test_reopening_gives_back_only_the_prediction_that_replaced_the_other(results_root):
    """
    O mesmo na etapa 2: `_manual_upload` troca a lista do método, não acrescenta a ela.

    Aqui o arquivo abandonado não muda a espécie: ele volta como predição, e o POA1 conta os
    epítopos dele junto com os do arquivo que a pessoa escolheu em seu lugar.
    """
    ctx = backend.create_run(results_root, "Predicao trocada", date(2026, 10, 16))
    ref = _reference(ctx)
    old = backend.save_upload(ctx.inputs_dir, "others_0_velho.fasta",
                              b">NS1_DENV1_OTHER_1_5\nMKTAY\n", ctx=ctx, role="prediction",
                              method="x", original_name="velho.fasta")
    new = backend.save_upload(ctx.inputs_dir, "others_0_novo.fasta",
                              b">NS1_DENV1_OTHER_3_7\nTAYIA\n", ctx=ctx, role="prediction",
                              method="x", original_name="novo.fasta")

    at = _start(run_path=str(ctx.root), ref_raw=[ref], poa1_sources={"x": [old]})
    at = at.sidebar.radio[0].set_value(STEP2).run()
    at.session_state["poa1_sources"] = {"x": [new]}
    at.run()
    assert not at.exception

    at = _reopen()

    assert not at.exception
    assert [Path(p).name for p in at.session_state["poa1_sources"]["x"]] == ["others_0_novo.fasta"]


def test_switching_analysis_does_not_carry_over_the_previous_forms(results_root):
    """
    Trocar de análise mantém os filtros de propósito, mas o que a nova análise gravou vence.

    `_reset_analysis_state` guarda os parâmetros entre análises porque abrir uma segunda costuma
    ser rodar os mesmos valores sobre outros dados. O que não pode é a segunda mostrar os alelos
    da primeira ao lado dos HTML dela própria, que saíram de outros.
    """
    first = backend.create_run(results_root, "Alelos A", date(2026, 10, 16))
    second = backend.create_run(results_root, "Alelos B", date(2026, 10, 17))   # a mais recente
    _step_two_with_alleles(second, "HLA-DPA1*01:03/DPB1*02:01", 21)

    at = _step_two_with_alleles(first, "HLA-DRB1*15:01", 12)
    at = next(b for b in at.sidebar.button if "Trocar" in b.label).click().run()
    at = next(b for b in at.button if "Abrir" in b.label).click().run()
    at = at.sidebar.radio[0].set_value(STEP2).run()

    assert not at.exception
    assert at.session_state["run_path"] == str(second.root)
    assert at.session_state["params_mhcii"] == {"alleles": "HLA-DPA1*01:03/DPB1*02:01",
                                                "length": 21}
    assert _widget(at.text_input, "Alelos").value == "HLA-DPA1*01:03/DPB1*02:01"
