"""Streamlit GUI for the POA pipeline — guided 6-step flow.

Run with:  streamlit run poa/gui/app.py   (or:  python run_gui.py)
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

# When launched via `streamlit run poa/gui/app.py`, Streamlit puts this file's directory on
# sys.path (not the project root), so `import poa` would fail. Add the repo root explicitly.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO

from poa.gui import backend, viz
from poa.logging_conf import setup_logging
from poa.services import (
    antigenic_client,
    bepipred_client,
    conservancy_client,
    mhcii_client,
    tmhmm_client,
)
from poa.services.base import ServiceError, ServiceUnavailable

setup_logging()
st.set_page_config(page_title="POA — Pipeline de Otimização de Antígenos", page_icon="🧬", layout="wide")


# --------------------------------------------------------------------------- session state
def _init_state():
    ss = st.session_state
    if "ctx" not in ss:
        ss.ctx = backend.WorkContext(Path(tempfile.mkdtemp(prefix="poa_gui_")))
    ss.setdefault("proteins_path", None)
    ss.setdefault("params_poa1", {})
    ss.setdefault("poa1_files", {})       # keys b2/b3/p/n/m/x -> path
    ss.setdefault("poa1_result", None)
    ss.setdefault("threshold", 70)
    ss.setdefault("conservancy_ready", False)
    ss.setdefault("poa2_result", None)
    return ss


ss = _init_state()
ctx: backend.WorkContext = ss.ctx


def _read_proteins():
    if not ss.proteins_path:
        return []
    return [(rec.id, str(rec.seq)) for rec in SeqIO.parse(ss.proteins_path, "fasta")]


# --------------------------------------------------------------------------- sidebar
st.sidebar.title("🧬 POA")
st.sidebar.caption("Pipeline de Otimização de Antígenos")

STEPS = [
    "1 · Proteínas & Parâmetros",
    "2 · Predições (auto/manual)",
    "3 · POA1 & Relatório",
    "4 · Conservancy Analysis",
    "5 · POA2 (conservação + topologia)",
    "6 · Resultados",
    "7 · Visualização (2D/3D)",
]


def _status_icon(done: bool) -> str:
    return "✅" if done else "⬜"


with st.sidebar:
    st.markdown("### Progresso")
    st.write(f"{_status_icon(bool(ss.proteins_path))} Proteínas carregadas")
    st.write(f"{_status_icon(bool(ss.poa1_files))} Resultados de predição")
    st.write(f"{_status_icon(ss.poa1_result is not None)} POA1 executado")
    st.write(f"{_status_icon(ss.conservancy_ready)} Conservancy pronta")
    st.write(f"{_status_icon(ss.poa2_result is not None)} POA2 executado")
    st.divider()
    step = st.radio("Etapa", STEPS)
    st.caption(f"Diretório de trabalho:\n`{ctx.root}`")


st.title("Pipeline de Otimização de Antígenos (POA)")


# =========================================================================== STEP 1
def step_inputs():
    st.header("1 · Proteínas e parâmetros do POA1")
    st.write("Envie o arquivo FASTA com as proteínas/poliproteínas usadas nas predições "
             "(cabeçalho no formato `Proteína_Espécie_ID`).")
    up = st.file_uploader("FASTA de proteínas (-f)", type=["fasta", "fa", "faa", "txt"])
    if up is not None:
        path = backend.save_upload(ctx.inputs_dir, "proteins.fasta", up.getvalue())
        ss.proteins_path = path
        st.success(f"Proteínas carregadas: {up.name}")

    if ss.proteins_path:
        seqs = _read_proteins()
        st.info(f"{len(seqs)} sequência(s) detectada(s): " + ", ".join(s[0] for s in seqs[:10]) +
                (" …" if len(seqs) > 10 else ""))

    st.subheader("Parâmetros de filtragem (opcionais)")
    c1, c2, c3 = st.columns(3)
    with c1:
        bmin = st.number_input("Bepipred min (bmin)", 0, 100, ss.params_poa1.get("bmin", 0))
        bmax = st.number_input("Bepipred max (bmax)", 0, 100, ss.params_poa1.get("bmax", 0))
    with c2:
        pmin = st.number_input("PAP min (pmin)", 0, 100, ss.params_poa1.get("pmin", 0))
        pmax = st.number_input("PAP max (pmax)", 0, 100, ss.params_poa1.get("pmax", 0))
    with c3:
        xmin = st.number_input("Outros min (xmin)", 0, 100, ss.params_poa1.get("xmin", 0))
        xmax = st.number_input("Outros max (xmax)", 0, 100, ss.params_poa1.get("xmax", 0))
    c4, c5, c6 = st.columns(3)
    with c4:
        mhla = st.selectbox("Classe HLA (MHC-II)", ["DR", "DP", "DQ"],
                            index=["DR", "DP", "DQ"].index(ss.params_poa1.get("mhla", "DR")))
    with c5:
        mic = st.number_input("IC50 threshold (mic, nM)", 0, 50000, ss.params_poa1.get("mic", 50))
    with c6:
        export = st.checkbox("Exportar .xlsx por método (-e)", value=ss.params_poa1.get("e", "n") == "y")

    ss.params_poa1 = dict(bmin=bmin, bmax=bmax, pmin=pmin, pmax=pmax, xmin=xmin, xmax=xmax,
                          mhla=mhla, mic=mic, e="y" if export else "n")
    st.success("Parâmetros salvos.") if ss.proteins_path else st.warning("Envie o FASTA de proteínas para continuar.")


# =========================================================================== STEP 2
def _automate_bepipred3():
    with st.status("Executando BepiPred-3.0 (local)…", expanded=True) as status:
        try:
            res = bepipred_client.predict(ss.proteins_path, cache=ctx.cache())
            path = backend.save_upload(ctx.inputs_dir, "bepipred3.fasta", res.content.encode("utf-8"))
            ss.poa1_files["b3"] = path
            ss.poa1_files.pop("b2", None)
            status.update(label=f"BepiPred-3.0 concluído ({res.source}).", state="complete")
        except (ServiceUnavailable, ServiceError) as exc:
            status.update(label="BepiPred-3.0 indisponível.", state="error")
            st.warning(f"{exc}\n\nUse o upload manual abaixo (plano B).")


def _automate_antigenic():
    with st.status("Executando EMBOSS antigenic (local)…", expanded=True) as status:
        try:
            res = antigenic_client.predict(ss.proteins_path, cache=ctx.cache())
            pap_txt = antigenic_client.to_pap_txt(res.content, ss.proteins_path)
            path = backend.save_upload(ctx.inputs_dir, "pap_antigenic.txt", pap_txt.encode("utf-8"))
            ss.poa1_files["p"] = path
            status.update(label=f"EMBOSS antigenic concluído ({res.source}).", state="complete")
        except (ServiceUnavailable, ServiceError) as exc:
            status.update(label="EMBOSS antigenic indisponível.", state="error")
            st.warning(f"{exc}\n\nUse o upload manual abaixo (plano B).")


def _automate_mhcii(alleles: str, length: int):
    seqs = _read_proteins()
    if not seqs:
        st.warning("Carregue o FASTA de proteínas primeiro.")
        return
    ok, fail = 0, 0
    with st.status("Consultando a API do IEDB (MHC-II) por proteína…", expanded=True) as status:
        for rec_id, seq in seqs:
            parts = rec_id.upper().split("_")
            protein = parts[0] if parts else rec_id.upper()
            specie = parts[1] if len(parts) > 1 else "SP"
            try:
                res = mhcii_client.predict(seq, alleles, length=length, cache=ctx.cache())
                mhcii_client.write_result_file(res.content, protein, specie, str(ctx.mhcii_dir))
                st.write(f"• {protein}_{specie}: ok ({res.source})")
                ok += 1
            except (ServiceUnavailable, ServiceError) as exc:
                st.write(f"• {protein}_{specie}: falhou — {exc}")
                fail += 1
        if ok:
            ss.poa1_files["m"] = str(ctx.mhcii_dir)
        state = "complete" if fail == 0 else ("error" if ok == 0 else "running")
        status.update(label=f"MHC-II: {ok} ok, {fail} falha(s).", state=state)


def _import_bepipred2_original():
    """Adapter for original manual-workflow BepiPred-2.0 JSONs (generic 'Sequence' antigen key)."""
    st.markdown("---")
    st.markdown("**Dados originais do fluxo manual (BepiPred-2.0)**")
    st.caption("Para JSONs do BepiPred-2.0 com chave de antígeno genérica (`Sequence`): reescrevemos "
               "o cabeçalho para `Proteína_Espécie_ID` e reconstruímos a proteína de referência (-f) "
               "a partir do próprio array `AA` do JSON.")
    c1, c2 = st.columns(2)
    with c1:
        specie = st.text_input("Espécie (ex.: DENV1, CHIKV)", "", key="rd_specie")
    with c2:
        protein = st.text_input("Proteína (ex.: E, E1)", "E", key="rd_protein")
    raw = st.file_uploader("BepiPred-2.0 original (.json)", type=["json"], key="up_b2_raw")
    if raw is not None and st.button("Adaptar e usar", key="rd_import_btn"):
        if not specie.strip():
            st.warning("Informe a espécie antes de adaptar.")
            return
        prep = backend.import_bepipred2(ctx, raw.name, raw.getvalue(),
                                        specie.strip(), (protein.strip() or "E"))
        ss.poa1_files["b2"] = prep.b2_json
        ss.poa1_files.pop("b3", None)
        adopted_ref = False
        if not ss.proteins_path:
            ss.proteins_path = prep.reference_fasta
            adopted_ref = True
        st.success(f"JSON adaptado como `{prep.header}` (referência de {len(prep.sequence)} aa).")
        if adopted_ref:
            st.info("Proteína de referência reconstruída do JSON e definida como -f (etapa 1).")
        if prep.has_x:
            st.warning("A referência contém o resíduo ambíguo 'X' — o POA1 vai recusá-la "
                       "(checagem de integridade do -f) até que ele seja resolvido.")


def _manual_upload(label: str, key: str, types, filename: str):
    up = st.file_uploader(label, type=types, key=f"up_{key}")
    if up is not None:
        path = backend.save_upload(ctx.inputs_dir, filename, up.getvalue())
        ss.poa1_files[key] = path
        st.success(f"{label}: {up.name} carregado.")


def step_predictions():
    st.header("2 · Predições de epítopos")
    if not ss.proteins_path:
        st.warning("Volte à etapa 1 e carregue o FASTA de proteínas.")
        return
    st.caption("Automatize quando houver via limpa; caso o serviço falhe, use o upload manual (plano B). "
               "Pelo menos um método é obrigatório.")

    # B cell — BepiPred
    with st.expander("🅱️ Células B — BepiPred", expanded=True):
        cols = st.columns([1, 2])
        with cols[0]:
            if st.button("Automatizar BepiPred-3.0 (local)"):
                _automate_bepipred3()
        with cols[1]:
            st.caption("BepiPred-3.0 via pacote local `bp3`. BepiPred-2.0: use upload do JSON.")
        _manual_upload("Upload BepiPred-2.0 (.json)", "b2", ["json"], "bepipred2.json")
        _manual_upload("Upload BepiPred-3.0 (.fasta)", "b3", ["fasta", "fa", "txt"], "bepipred3.fasta")
        _import_bepipred2_original()

    # B cell — PAP/IMED via EMBOSS
    with st.expander("🅱️ Células B — PAP/IMED (antigenicidade)", expanded=False):
        if st.button("Automatizar via EMBOSS antigenic (local)"):
            _automate_antigenic()
        _manual_upload("Upload PAP/IMED (.txt)", "p", ["txt"], "pap.txt")

    # Tc — NetCTL
    with st.expander("🅃 Células T-citotóxicas — NetCTL 1.2", expanded=False):
        st.caption("Sem API/binário Windows. Automação por navegador é opcional (etapa futura). "
                   "Plano B: faça a predição no site e envie o .html.")
        _manual_upload("Upload NetCTL (.html)", "n", ["html", "htm"], "netctl.html")

    # Th — MHC-II
    with st.expander("🅃 Células T-auxiliares — MHC-II Binding (IEDB)", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            alleles = st.text_input("Alelos (separados por vírgula)", "HLA-DRB1*01:01")
        with c2:
            length = st.number_input("Comprimento do peptídeo", 9, 25, 15)
        if st.button("Automatizar MHC-II (API IEDB)"):
            _automate_mhcii(alleles, int(length))
        st.caption("Plano B: envie os .html (um por proteína, nomeados `Proteína_Espécie.html`).")
        ups = st.file_uploader("Upload MHC-II (.html, múltiplos)", type=["html", "htm"],
                               accept_multiple_files=True, key="up_m_multi")
        if ups:
            for up in ups:
                backend.save_upload(ctx.mhcii_dir, up.name, up.getvalue())
            ss.poa1_files["m"] = str(ctx.mhcii_dir)
            st.success(f"{len(ups)} arquivo(s) MHC-II carregado(s).")

    # Others
    with st.expander("➕ Outros preditores (FASTA padronizado)", expanded=False):
        _manual_upload("Upload outros (.fasta)", "x", ["fasta", "fa", "txt"], "others.fasta")

    st.divider()
    if ss.poa1_files:
        st.success("Métodos prontos: " + ", ".join(sorted(ss.poa1_files.keys())))
    else:
        st.warning("Nenhum resultado de predição ainda.")


# =========================================================================== STEP 3
def step_poa1():
    st.header("3 · POA1 — consolidação, ranqueamento e relatório")
    if not ss.proteins_path or not ss.poa1_files:
        st.warning("Complete as etapas 1 e 2 (proteínas + ao menos um método).")
        return
    if st.button("▶️ Executar POA1", type="primary"):
        with st.status("Executando POA1…", expanded=True) as status:
            try:
                result = backend.run_poa1(ctx, ss.poa1_files, ss.params_poa1, ss.proteins_path)
                ss.poa1_result = result
                status.update(label="POA1 concluído.", state="complete")
            except Exception as exc:
                status.update(label="POA1 falhou.", state="error")
                st.exception(exc)
                return

    if ss.poa1_result is not None:
        df = ss.poa1_result.predictions
        st.subheader("Epítopos consolidados")
        st.dataframe(df, use_container_width=True, height=320)

        summ = backend.predictions_summary(df)
        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(px.bar(summ["by_method"], x="Method", y="count",
                                   title="Epítopos por método"), use_container_width=True)
        with c2:
            st.plotly_chart(px.bar(summ["by_species"], x="Specie", y="count",
                                   title="Epítopos por espécie"), use_container_width=True)
        st.plotly_chart(px.bar(summ["length_hist"], x="length", y="count",
                               title="Distribuição de comprimento dos epítopos"), use_container_width=True)

        report = Path(ss.poa1_result.report_path)
        if report.exists():
            with st.expander("📄 Relatório de análise (Analysis_report.txt)"):
                st.code(report.read_text(encoding="utf-8"))


# =========================================================================== STEP 4
def step_conservancy():
    st.header("4 · Epitope Conservancy Analysis")
    if ss.poa1_result is None:
        st.warning("Execute o POA1 (etapa 3) primeiro.")
        return
    st.caption("O IEDB não oferece API nem versão local desta ferramenta. O POA reimplementa o "
               "cálculo localmente (Bui et al., 2007). Você também pode enviar os CSVs manualmente.")

    ss.threshold = st.number_input("Limiar de identidade de sequência (%)", 1, 100, ss.threshold)
    mode = st.radio("Modo", ["Reimplementação local (recomendado)", "Upload manual dos CSVs"])

    if mode.startswith("Reimplementação"):
        st.write("Conjunto de proteínas para comparação (padrão: as proteínas carregadas na etapa 1).")
        comp = st.file_uploader("FASTA de comparação (opcional)", type=["fasta", "fa", "txt"], key="comp_fasta")
        comp_path = ss.proteins_path
        if comp is not None:
            comp_path = backend.save_upload(ctx.inputs_dir, "comparison.fasta", comp.getvalue())
        if st.button("▶️ Calcular conservância (local)", type="primary"):
            with st.status("Calculando conservância…", expanded=True) as status:
                try:
                    written = conservancy_client.run_conservancy_for_dir(
                        str(ctx.conservancy_epitopes_dir), comp_path, float(ss.threshold),
                        str(ctx.conservancy_csv_dir))
                    ss.conservancy_ready = len(written) > 0
                    status.update(label=f"{len(written)} CSV(s) gerado(s).", state="complete")
                except Exception as exc:
                    status.update(label="Falha no cálculo de conservância.", state="error")
                    st.exception(exc)
    else:
        ups = st.file_uploader("CSVs da Conservancy Analysis (IEDB)", type=["csv"],
                               accept_multiple_files=True, key="cons_csvs")
        if ups:
            for up in ups:
                backend.save_upload(ctx.conservancy_csv_dir, up.name, up.getvalue())
            ss.conservancy_ready = True
            st.success(f"{len(ups)} CSV(s) carregado(s).")

    csvs = list(ctx.conservancy_csv_dir.glob("*.csv"))
    if csvs:
        st.info("CSVs disponíveis: " + ", ".join(c.name for c in csvs))


# =========================================================================== STEP 5
def step_poa2():
    st.header("5 · POA2 — conservação e topologia de membrana")
    if not ss.conservancy_ready:
        st.warning("Prepare a Conservancy Analysis (etapa 4) primeiro.")
        return
    if not tmhmm_client.is_available():
        st.error("pyTMHMM não está instalado. Instale com `pip install pyTMHMM` para rodar o POA2.")
        return

    c1, c2, c3 = st.columns(3)
    with c1:
        objective = st.radio("Objetivo", ["conserved", "unique"],
                             format_func=lambda x: "Conservados (≥)" if x == "conserved" else "Únicos (<)")
        rf = st.selectbox("FASTA por topologia (-rf)",
                          [None, 0, 1, 2, 3],
                          format_func=lambda v: {None: "não gerar", 0: "todos", 1: "externos",
                                                 2: "transmembrana", 3: "internos"}[v])
    with c2:
        imin = st.number_input("Identidade mín. (%)", 0, 100, 60)
        imax = st.number_input("Identidade máx. (%)", 0, 100, 100)
    with c3:
        mmatch = st.number_input("% de sequências com match (-m)", 0, 100, 60)
        st.metric("Limiar (t)", ss.threshold)

    if st.button("▶️ Executar POA2", type="primary"):
        params = dict(objective=objective, t=ss.threshold, imin=imin, imax=imax, m=mmatch, rf=rf)
        with st.status("Executando POA2 (conservância + TMHMM)…", expanded=True) as status:
            try:
                ss.poa2_result = backend.run_poa2(ctx, params, ss.proteins_path)
                status.update(label="POA2 concluído.", state="complete")
            except Exception as exc:
                status.update(label="POA2 falhou.", state="error")
                st.exception(exc)


# =========================================================================== STEP 6
def _download(path: Path, label: str, mime: str):
    if path and path.exists():
        st.download_button(label, data=path.read_bytes(), file_name=path.name, mime=mime)


def step_results():
    st.header("6 · Resultados finais")
    if ss.poa2_result is None:
        st.warning("Execute o POA2 (etapa 5) para ver os resultados finais.")
        # still offer POA1 downloads
    if ss.poa2_result is not None:
        df = ss.poa2_result.results
        st.subheader("Epítopos selecionados (POA2)")
        st.dataframe(df, use_container_width=True, height=320)

        topo = backend.topology_summary(df)
        if not topo.empty:
            st.plotly_chart(px.bar(topo, x="region", y="mean_fraction",
                                   title="Fração média por topologia de membrana"),
                            use_container_width=True)
        st.subheader("Downloads")
        _download(Path(ss.poa2_result.xlsx_path), "⬇️ Planilha POA2 (.xlsx)",
                  "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        if ss.poa2_result.fasta_path:
            _download(Path(ss.poa2_result.fasta_path), "⬇️ FASTA POA2 (.fasta)", "text/plain")

    if ss.poa1_result is not None:
        st.subheader("Saídas do POA1")
        _download(Path(ss.poa1_result.report_path), "⬇️ Relatório (.txt)", "text/plain")
        for fasta in sorted(ctx.conservancy_epitopes_dir.glob("*.fasta")):
            _download(fasta, f"⬇️ {fasta.name}", "text/plain")


# =========================================================================== STEP 7
def step_viz():
    st.header("7 · Visualização dos epítopos (2D / 3D)")
    if ss.poa1_result is None:
        st.warning("Execute o POA1 (etapa 3) para ter epítopos a visualizar.")
        return
    df = ss.poa1_result.predictions

    # --- 2D epitope map (no structure required) ---
    st.subheader("Mapa 2D de epítopos")
    st.caption("Cada segmento é um epítopo ao longo da sequência; uma trilha por proteína/espécie, "
               "colorido por método.")
    st.plotly_chart(viz.epitope_map_figure(df), use_container_width=True)

    # --- 3D structure viewer (upload a PDB) ---
    st.subheader("Estrutura 3D (envie um PDB)")
    st.caption("Estilo Discovery Studio: epítopos destacados sobre a estrutura. O visualizador usa "
               "3Dmol.js (precisa de internet para renderizar).")

    tracks = sorted({f"{r.Protein}|{r.Specie}" for r in df.itertuples()})
    choice = st.selectbox("Epítopos de qual proteína/espécie destacar?", tracks)
    c1, c2, c3 = st.columns(3)
    with c1:
        base_style = st.selectbox("Representação", ["cartoon", "stick", "sphere"], index=0)
    with c2:
        surface = st.checkbox("Mostrar superfície", value=False)
    with c3:
        chain = st.text_input("Cadeia (opcional)", "")

    pdb_up = st.file_uploader("Arquivo .pdb", type=["pdb", "ent"], key="pdb_up")
    if pdb_up is not None:
        prot, spec = choice.split("|")
        ranges = viz.epitope_ranges(df, protein=prot, specie=spec)
        try:
            html = viz.build_3dmol_view_html(
                pdb_up.getvalue().decode("utf-8", "replace"),
                ranges,
                chain=chain or None,
                base_style=base_style,
                show_surface=surface,
            )
            components.html(html, height=520)
            st.caption(f"{len(ranges)} epítopo(s) de {prot}_{spec} destacado(s) em vermelho. "
                       "Atenção: a numeração de resíduo do PDB precisa corresponder à posição na "
                       "sequência usada nas predições (cuidado com gaps/offset).")
        except Exception as exc:  # noqa: BLE001 - surface the error to the user, don't crash the app
            st.error(f"Não foi possível renderizar a estrutura: {exc}")
    else:
        st.info("Envie um arquivo .pdb para ver os epítopos destacados na estrutura 3D.")


# --------------------------------------------------------------------------- dispatch
_DISPATCH = {
    STEPS[0]: step_inputs,
    STEPS[1]: step_predictions,
    STEPS[2]: step_poa1,
    STEPS[3]: step_conservancy,
    STEPS[4]: step_poa2,
    STEPS[5]: step_results,
    STEPS[6]: step_viz,
}
_DISPATCH[step]()
