"""Streamlit GUI for the POA pipeline — guided 6-step flow.

Run with:  streamlit run poa/gui/app.py   (or:  python run_gui.py)
"""
from __future__ import annotations

import sys
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

from poa.core.parsers import conservancy as conservancy_parser
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
#: session keys holding data that belongs to one analysis — cleared when switching analyses,
#: so a new run never inherits the previous one's file paths or results
_ANALYSIS_KEYS = ("proteins_path", "ref_raw", "ref_uploads", "poa1_sources", "prepared_b2",
                  "poa1_result", "conservancy_ready", "comparison_sets", "poa2_result")


def _init_state():
    ss = st.session_state
    ss.setdefault("run_path", None)        # the analysis folder this session writes to
    ss.setdefault("proteins_path", None)   # merged -f reference (all species)
    ss.setdefault("ref_raw", [])           # protein FASTAs as uploaded in step 1
    ss.setdefault("ref_swap", False)       # source headers are Specie_Protein -> swap them
    ss.setdefault("ref_uploads", [])       # the same files after optional header normalisation
    ss.setdefault("params_poa1", {})
    # one *list* of files per method (b2/b3/p/n/x) so several species can be analysed together;
    # 'm' (MHC-II) stays a single directory of Protein_Specie.html files
    ss.setdefault("poa1_sources", {})
    ss.setdefault("prepared_b2", [])       # Prepared records from the BepiPred-2.0 adapter
    ss.setdefault("poa1_result", None)
    ss.setdefault("threshold", 70)
    ss.setdefault("cons_operator", ">=")   # '>=' conserved | '<' unique — must match the POA2 objective
    ss.setdefault("conservancy_ready", False)
    ss.setdefault("comparison_sets", {})   # SPECIE -> comparison FASTA path
    ss.setdefault("poa2_result", None)
    return ss


ss = _init_state()
ctx: backend.WorkContext = None   # bound by _select_run() before any step renders


def _collect_sources():
    """All per-species prediction files, keyed by POA1 method argument."""
    src = {k: list(v) for k, v in ss.poa1_sources.items() if v}
    if ss.prepared_b2:
        src.setdefault("b2", []).extend(pr.b2_json for pr in ss.prepared_b2)
    return src


def _refresh_reference(allow_clear: bool = False):
    """
    Merge every reference FASTA (uploaded + rebuilt by the adapter) into a single ``-f``.

    With nothing to merge the current ``-f`` is kept unless ``allow_clear`` — step 1 re-renders on
    every interaction, and a reference adopted elsewhere (e.g. resumed from ``results/`` in step 5)
    must not be wiped just because no file is sitting in the uploader.
    """
    parts = list(ss.ref_uploads) + [pr.reference_fasta for pr in ss.prepared_b2]
    if not parts:
        if allow_clear:
            ss.proteins_path = None
        return ss.proteins_path
    if len(parts) == 1:
        ss.proteins_path = parts[0]
    else:
        ss.proteins_path = backend.merge_fastas(parts, str(ctx.inputs_dir / "proteins_all.fasta"))
    return ss.proteins_path


def _read_proteins():
    if not ss.proteins_path:
        return []
    return [(rec.id, str(rec.seq)) for rec in SeqIO.parse(ss.proteins_path, "fasta")]


# --------------------------------------------------------------------------- analysis folder
def _reset_analysis_state():
    """Drop everything tied to the previous analysis — its paths point into another folder."""
    for key in _ANALYSIS_KEYS:
        del ss[key]
    _init_state()


def _adopt_run(path: Path):
    ss.run_path = str(path)
    _reset_analysis_state()
    st.rerun()


def _select_run() -> backend.WorkContext:
    """
    Bind this session to one analysis folder, asking which before anything else runs.

    Every analysis gets its own ``results/DDMMAAAA-NOME_DO_TESTE``, so a second run of the same
    data never overwrites the first and the two can be compared side by side.
    """
    root = backend.default_results_root()

    if ss.run_path:
        path = Path(ss.run_path)
        with st.sidebar:
            st.markdown("### Análise")
            st.success(f"**{path.name or root.name}**")
            if st.button("Trocar de análise"):
                ss.run_path = None
                st.rerun()
        return backend.open_run(root, path)

    runs = backend.list_runs(root)
    st.title("Pipeline de Otimização de Antígenos (POA)")
    st.header("Escolha a análise")
    st.caption(f"Cada análise fica na sua própria pasta dentro de `{root}`, nomeada "
               "`DDMMAAAA-NOME_DO_TESTE`. Assim uma nova execução não sobrescreve a anterior e dá "
               "para comparar as duas.")

    tab_new, tab_open = st.tabs(["Nova análise", f"Abrir existente ({len(runs)})"])
    with tab_new:
        label = st.text_input("Nome do teste", placeholder="ex.: DENV world set")
        if label.strip():
            st.caption(f"Pasta: `{backend.run_dir_name(label)}`")
        if st.button("Criar análise", type="primary", disabled=not label.strip()):
            _adopt_run(backend.create_run(root, label).root)
    with tab_open:
        if not runs:
            st.info("Nenhuma análise encontrada ainda. Crie a primeira na outra aba.")
        else:
            def _fmt(i):
                run = runs[i]
                when = run.date.strftime("%d/%m/%Y") if run.date else "sem data"
                return f"{run.label} — {when} — {run.progress}"

            picked = st.selectbox("Análises gravadas", range(len(runs)), format_func=_fmt)
            chosen = runs[picked]
            st.code(str(chosen.path), language=None)
            if chosen.legacy:
                st.info("Esta é a pasta de resultados anterior à convenção de pastas datadas. Ela "
                        "continua utilizável onde está — nada foi movido.")
            if st.button("Abrir análise"):
                _adopt_run(chosen.path)
    st.stop()


# --------------------------------------------------------------------------- sidebar
st.sidebar.title("🧬 POA")
st.sidebar.caption("Pipeline de Otimização de Antígenos")

ctx = _select_run()

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
    st.write(f"{_status_icon(bool(_collect_sources()))} Resultados de predição")
    st.write(f"{_status_icon(ss.poa1_result is not None)} POA1 executado")
    st.write(f"{_status_icon(ss.conservancy_ready)} Conservancy pronta")
    st.write(f"{_status_icon(ss.poa2_result is not None)} POA2 executado")
    st.divider()
    step = st.radio("Etapa", STEPS)
    st.divider()
    st.markdown("### Resultados")
    st.caption("Tudo desta análise é gravado aqui. Outra análise, outra pasta — repetir uma "
               "execução nesta sobrescreve os arquivos de mesmo nome:")
    st.code(str(ctx.root), language=None)


st.title("Pipeline de Otimização de Antígenos (POA)")


# =========================================================================== STEP 1
def step_inputs():
    st.header("1 · Proteínas e parâmetros do POA1")
    st.write("Envie o(s) arquivo(s) FASTA com as proteínas/poliproteínas usadas nas predições "
             "(cabeçalho no formato `Proteína_Espécie_ID`, ex.: `E_DENV1_ref`). "
             "**Para analisar mais de uma espécie**, inclua todas as proteínas — no mesmo FASTA "
             "ou em vários arquivos; eles são unidos em um único `-f`.")
    ups = st.file_uploader("FASTA de proteínas (-f)", type=["fasta", "fa", "faa", "txt"],
                           accept_multiple_files=True)
    if ups:
        ss.ref_raw = [backend.save_upload(ctx.inputs_dir, f"ref_{i}_{u.name}", u.getvalue())
                      for i, u in enumerate(ups)]

    swap = st.checkbox("Meus cabeçalhos estão em `Espécie_Proteína` (ex.: `denv1_ns1`) — inverter "
                       "os dois primeiros campos", value=ss.ref_swap, key="ref_swap_cb")
    ss.ref_swap = swap

    if ss.ref_raw:
        if swap:
            ss.ref_uploads = [backend.normalize_fasta_headers(
                ss.ref_raw, str(ctx.inputs_dir / "proteins_normalized.fasta"), swap=True)]
        else:
            ss.ref_uploads = list(ss.ref_raw)
        st.success(f"{len(ss.ref_raw)} arquivo(s) de proteínas carregado(s).")
    _refresh_reference()

    if ss.proteins_path:
        rows = backend.parse_fasta_headers(ss.proteins_path)
        species = backend.species_in_fasta(ss.proteins_path)
        total_aa = sum(r["length"] for r in rows)
        st.write(f"**{len(rows)} sequência(s)**, {total_aa} resíduos no total — veja como cada "
                 "cabeçalho foi interpretado e o que cada registro contém:")
        st.dataframe(pd.DataFrame([{"Cabeçalho": r["header"], "Proteína": r["protein"],
                                    "Espécie": r["specie"], "ID": r["id"],
                                    "Resíduos": r["length"],
                                    "Ambíguos": r["ambiguous"] or "—"} for r in rows]),
                     use_container_width=True, hide_index=True)

        empty = [r["header"] for r in rows if r["length"] == 0]
        if empty:
            st.error("Registro(s) sem sequência: " + ", ".join(empty) +
                     ". Um cabeçalho sem resíduos não serve de referência para a conservância.")
        with_x = [r["header"] for r in rows if "X" in r["ambiguous"]]
        if with_x:
            st.warning("Sequência(s) com o resíduo ambíguo `X`: " + ", ".join(with_x) +
                       ". O POA1 recusa epítopos que contenham `X` (checagem de integridade do "
                       "`-f`) — resolva antes de executar a etapa 3.")
        other_amb = sorted({c for r in rows for c in r["ambiguous"] if c != "X"})
        if other_amb:
            st.caption("Outros resíduos não-padrão presentes: " + ", ".join(other_amb) +
                       " — não bloqueiam o POA1, mas afetam o cálculo de identidade.")

        if not all(r["ok"] for r in rows):
            st.error("Há cabeçalhos sem separador `_`, dos quais nenhuma espécie pode ser extraída. "
                     "Use `Proteína_Espécie_ID` (ex.: `NS1_DENV1_ref`).")
        elif species:
            st.success(f"{len(species)} espécie(s) identificada(s): " + ", ".join(species))
            st.caption("Se a coluna Espécie mostrar o nome da proteína, marque a caixa de inversão acima.")
        else:
            st.error("Nenhuma espécie identificada nos cabeçalhos. Use o formato "
                     "`Proteína_Espécie_ID` (ex.: `NS1_DENV1_ref`) — sem ele o POA não consegue "
                     "separar as espécies nas etapas 3 e 4.")

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
    if ss.proteins_path:
        st.success("Parâmetros salvos.")
    else:
        st.warning("Envie o FASTA de proteínas para continuar.")


# =========================================================================== STEP 2
def _automate_bepipred3():
    with st.status("Executando BepiPred-3.0 (local)…", expanded=True) as status:
        try:
            res = bepipred_client.predict(ss.proteins_path, cache=ctx.cache())
            path = backend.save_upload(ctx.inputs_dir, "bepipred3.fasta", res.content.encode("utf-8"))
            ss.poa1_sources["b3"] = [path]
            ss.poa1_sources.pop("b2", None)
            ss.prepared_b2 = []
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
            ss.poa1_sources["p"] = [path]
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
            ss.poa1_sources["m"] = [str(ctx.mhcii_dir)]
        state = "complete" if fail == 0 else ("error" if ok == 0 else "running")
        status.update(label=f"MHC-II: {ok} ok, {fail} falha(s).", state=state)


def _import_bepipred2_original():
    """Adapter for original manual-workflow BepiPred-2.0 JSONs (non-conforming antigen keys).

    Handles both shapes of the real data: one JSON per species, and a single JSON whose antigens
    are several organisms (a submission of ``denv1_ns1`` + ``denv2_ns1`` comes back keyed
    ``denv1``/``denv2``, which the ``Protein_Specie_ID`` regex cannot resolve). Each antigen gets
    its own species/protein, and successive imports are accumulated.
    """
    st.markdown("**BepiPred-2.0 (.json)**")
    st.caption("O arquivo é verificado no envio. Se as chaves de antígeno já seguirem "
               "`Proteína_Espécie_ID`, ele é usado direto; caso contrário (`Sequence`, `denv1`, …) "
               "você mapeia cada antígeno abaixo — reescrevemos o cabeçalho de **cada um** e "
               "reconstruímos as proteínas de referência (-f) a partir do próprio array `AA` do "
               "JSON. Vários antígenos no mesmo arquivo podem ser espécies diferentes.")

    raw = st.file_uploader("BepiPred-2.0 (.json)", type=["json"], key="up_b2_raw")
    if raw is not None:
        raw_path = backend.save_upload(ctx.inputs_dir, f"raw_{raw.name}", raw.getvalue())
        try:
            keys = backend.antigen_keys(raw_path)
            bad = backend.nonconforming_antigen_keys(raw_path)
        except Exception as exc:  # noqa: BLE001 - malformed upload must not crash the app
            st.error(f"Não foi possível ler o JSON: {exc}")
            keys, bad = [], []

        if not keys:
            st.error("Nenhum bloco `antigens` encontrado no JSON.")
        elif not bad:
            # keys already conform -> the parser resolves species on its own, no adapter needed
            ss.poa1_sources["b2"] = [raw_path]
            ss.poa1_sources.pop("b3", None)
            st.success(f"{len(keys)} antígeno(s) com chave no padrão `Proteína_Espécie_ID`: " +
                       ", ".join(f"`{k}`" for k in keys[:6]) + (" …" if len(keys) > 6 else "") +
                       ". Arquivo usado diretamente.")
            st.info("Lembre-se de que o FASTA de proteínas (-f) da etapa 1 precisa conter essas "
                    "mesmas espécies.")
        else:
            ss.poa1_sources.pop("b2", None)   # the raw file is unusable as-is
            st.error("As chaves de antígeno deste JSON não seguem `Proteína_Espécie_ID`: " +
                     ", ".join(f"`{k}`" for k in bad[:6]) + (" …" if len(bad) > 6 else "") +
                     ". Sem mapeá-las, o POA1 não consegue identificar a espécie e todos os "
                     "epítopos ficam sem espécie. Preencha os campos abaixo.")
            st.write(f"**{len(keys)} antígeno(s) no arquivo.** Informe a espécie e a proteína de "
                     "cada um — o cabeçalho passa a ser `Proteína_Espécie_ID`.")
            mapping = {}
            for i, key in enumerate(keys):
                tokens = str(key).upper().split("_")
                vis = "visible" if i == 0 else "collapsed"
                c1, c2, c3 = st.columns([2, 2, 2])
                with c1:
                    st.text_input("Chave no JSON", str(key), disabled=True,
                                  key=f"rd_key_{i}", label_visibility=vis)
                with c2:
                    sp = st.text_input("Espécie", tokens[0], key=f"rd_sp_{i}", label_visibility=vis)
                with c3:
                    prot = st.text_input("Proteína", tokens[1] if len(tokens) > 1 else "E",
                                         key=f"rd_prot_{i}", label_visibility=vis)
                mapping[key] = ((prot.strip().upper() or "E"), sp.strip().upper(), "ref")

            if st.button("Adaptar e adicionar", key="rd_import_btn"):
                if any(not v[1] for v in mapping.values()):
                    st.warning("Informe a espécie de todos os antígenos antes de adaptar.")
                    return
                already = {sp for pr in ss.prepared_b2 for sp in pr.species}
                dup = sorted({v[1] for v in mapping.values()} & already)
                if dup:
                    st.warning(f"Já importado(s): {', '.join(dup)}. Limpe as importações antes de repetir.")
                    return
                first = mapping[keys[0]]
                prep = backend.import_bepipred2(ctx, raw.name, raw.getvalue(),
                                                specie=first[1], protein=first[0], mapping=mapping)
                ss.prepared_b2.append(prep)
                ss.poa1_sources.pop("b3", None)
                _refresh_reference()
                st.success("Antígenos adaptados: " +
                           ", ".join(f"`{e.header}` ({len(e.sequence)} aa)" for e in prep.entries))
                st.info("Proteínas de referência reconstruídas do JSON e somadas ao -f (etapa 1).")
                if prep.has_x:
                    st.warning("A referência contém o resíduo ambíguo 'X' — o POA1 vai recusá-la "
                               "(checagem de integridade do -f) até que ele seja resolvido.")

    if ss.prepared_b2:
        st.write("**Espécies importadas:** " +
                 ", ".join(f"`{e.header}`" for pr in ss.prepared_b2 for e in pr.entries))
        if st.button("Limpar espécies importadas", key="rd_clear_btn"):
            ss.prepared_b2 = []
            _refresh_reference(allow_clear=True)
            st.rerun()


def _manual_upload(label: str, key: str, types, filename: str):
    """Upload one file *per species* for a method; they are merged into the single POA1 input."""
    ups = st.file_uploader(label, type=types, key=f"up_{key}", accept_multiple_files=True)
    if ups:
        stem = filename.rsplit(".", 1)[0]
        ss.poa1_sources[key] = [
            backend.save_upload(ctx.inputs_dir, f"{stem}_{i}_{up.name}", up.getvalue())
            for i, up in enumerate(ups)
        ]
        st.success(f"{label}: {len(ups)} arquivo(s) — " + ", ".join(u.name for u in ups))


def _human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB"):
        if n < 1024 or unit == "MB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} MB"


def _render_source_previews(sources):
    """
    Show what each collected prediction file holds, before POA1 ever reads it.

    Step 2 is where a wrong file or an unreadable header convention is cheapest to fix; without
    this the first sign of either was an empty or mis-grouped table two steps later.
    """
    st.subheader("Conteúdo dos arquivos enviados")
    st.caption("Este é o **conteúdo bruto** de cada arquivo, lido do próprio arquivo. Os epítopos "
               "em si são extraídos pelo POA1 (etapa 3), com os parâmetros de filtragem da etapa 1 "
               "— por isso aqui aparece o que o arquivo traz, não o que o POA1 vai selecionar.")

    previews = backend.summarize_sources(sources)
    all_species = sorted({sp for p in previews for sp in p.species})
    if all_species:
        st.info("Espécies encontradas nos arquivos de predição: " + ", ".join(all_species))

    for prev in previews:
        unit = {"fasta": "registro(s)", "json": "antígeno(s)",
                "dir": "arquivo(s)", "text": "linha(s)"}[prev.kind]
        count = "?" if prev.n_records is None else prev.n_records
        head = f"{prev.label} · `{prev.name}` — {count} {unit} · {_human_bytes(prev.size_bytes)}"
        with st.expander(head, expanded=len(previews) == 1):
            if prev.species:
                st.write("**Espécies:** " + ", ".join(prev.species))
            if prev.proteins:
                st.write("**Proteínas:** " + ", ".join(prev.proteins))
            if prev.rows:
                st.dataframe(pd.DataFrame(prev.rows), use_container_width=True,
                             hide_index=True, height=240)
            if prev.head:
                st.code(prev.head, language=None)
            for note in prev.notes:
                st.warning(note)


def step_predictions():
    st.header("2 · Predições de epítopos")
    if not ss.proteins_path:
        st.warning("Volte à etapa 1 e carregue o FASTA de proteínas.")
        return
    st.caption("Automatize quando houver via limpa; caso o serviço falhe, use o upload manual (plano B). "
               "Pelo menos um método é obrigatório. **Multi-espécie:** envie um arquivo por espécie "
               "em cada método — eles são unidos automaticamente antes do POA1.")

    # B cell — BepiPred
    with st.expander("🅱️ Células B — BepiPred", expanded=True):
        cols = st.columns([1, 2])
        with cols[0]:
            if st.button("Automatizar BepiPred-3.0 (local)"):
                _automate_bepipred3()
        with cols[1]:
            st.caption("BepiPred-3.0 via pacote local `bp3`. BepiPred-2.0: use upload do JSON.")
        _manual_upload("Upload BepiPred-3.0 (.fasta)", "b3", ["fasta", "fa", "txt"], "bepipred3.fasta")
        st.markdown("---")
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
            ss.poa1_sources["m"] = [str(ctx.mhcii_dir)]
            st.success(f"{len(ups)} arquivo(s) MHC-II carregado(s).")

    # Others
    with st.expander("➕ Outros preditores (FASTA padronizado)", expanded=False):
        _manual_upload("Upload outros (.fasta)", "x", ["fasta", "fa", "txt"], "others.fasta")

    st.divider()
    sources = _collect_sources()
    if not sources:
        st.warning("Nenhum resultado de predição ainda.")
        return

    st.success("Métodos prontos: " +
               ", ".join(f"{k} ({len(v)} arquivo(s))" for k, v in sorted(sources.items())))
    _render_source_previews(sources)


# =========================================================================== STEP 3
def step_poa1():
    st.header("3 · POA1 — consolidação, ranqueamento e relatório")
    sources = _collect_sources()
    if not ss.proteins_path or not sources:
        st.warning("Complete as etapas 1 e 2 (proteínas + ao menos um método).")
        return

    ref_species = backend.species_in_fasta(ss.proteins_path)
    st.caption(f"Referência (-f): {len(ref_species)} espécie(s) — " +
               (", ".join(ref_species) or "nenhuma identificada"))
    if st.button("▶️ Executar POA1", type="primary"):
        with st.status("Executando POA1…", expanded=True) as status:
            try:
                files = backend.consolidate_sources(sources, ctx.inputs_dir)
                for key, path in sorted(files.items()):
                    n = len(sources[key])
                    st.write(f"• {key}: {Path(path).name}" + (f" ({n} arquivos unidos)" if n > 1 else ""))
                result = backend.run_poa1(ctx, files, ss.params_poa1, ss.proteins_path)
                ss.poa1_result = result
                ss.conservancy_ready = False
                status.update(label="POA1 concluído.", state="complete")
            except Exception as exc:
                status.update(label="POA1 falhou.", state="error")
                st.exception(exc)
                return

    if ss.poa1_result is not None:
        df = ss.poa1_result.predictions

        diag = backend.species_diagnostics(ss.proteins_path, df)
        if diag["predicted"]:
            st.success("Espécies identificadas nos epítopos: " + ", ".join(diag["predicted"]))
        else:
            st.error("Nenhuma espécie foi identificada nos epítopos.")
        if diag["unnamed"]:
            st.error("Há epítopos sem espécie no cabeçalho. Os arquivos de predição precisam usar "
                     "`Proteína_Espécie_ID`; para JSONs originais do BepiPred-2.0 use o adaptador "
                     "da etapa 2.")
        if diag["missing_from_reference"]:
            st.error("Espécies presentes nas predições mas **ausentes do FASTA de referência (-f)**: "
                     + ", ".join(diag["missing_from_reference"]) +
                     ". Adicione as proteínas dessas espécies na etapa 1 — sem elas a análise de "
                     "conservância (etapa 4) compara epítopos com as proteínas erradas.")
        if diag["unpredicted"]:
            st.warning("Espécies no -f sem nenhum epítopo predito: " + ", ".join(diag["unpredicted"]))
        if diag["reference"] and diag["predicted"] and not set(diag["reference"]) & set(diag["predicted"]):
            st.warning("Nenhuma espécie do -f coincide com as das predições. Verifique na etapa 1 a "
                       "tabela de leitura dos cabeçalhos — se a coluna Espécie estiver mostrando o "
                       "nome da proteína, marque a caixa de inversão `Espécie_Proteína`.")

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
def _declare_uploaded_threshold():
    """
    Record the threshold/operator behind manually uploaded conservancy CSVs.

    POA2 validates what it is filtering against this record. IEDB spells both values in the CSV
    header, so they are read from there when present and only asked for when they are not.
    """
    csv_dir = str(ctx.conservancy_csv_dir)
    source = conservancy_parser.inspect_directory(csv_dir)

    if source.origin == "meta":
        st.success(f"Registrado: limiar **{conservancy_parser.format_threshold(source.threshold)}%**, "
                   f"critério **{source.operator}**. A etapa 5 vai conferir estes valores.")
        return

    if source.origin == "header":
        detected_op = source.operator or ">="
        st.info(f"O cabeçalho dos CSVs declara limiar "
                f"**{conservancy_parser.format_threshold(source.threshold)}%** e critério "
                f"**{detected_op}**.")
        if st.button("Registrar estes valores", key="cons_adopt_header"):
            conservancy_parser.write_metadata(csv_dir, source.threshold, detected_op, source="iedb")
            ss.threshold = int(source.threshold)
            ss.cons_operator = ">=" if conservancy_parser.operator_family(detected_op) == "ge" else "<"
            st.rerun()
        return

    st.warning("O cabeçalho destes CSVs não traz o limiar usado. Informe-o abaixo — sem isso a "
               "etapa 5 não tem como validar o que está filtrando.")
    c1, c2 = st.columns([1, 2])
    with c1:
        man_t = st.number_input("Limiar usado (%)", 1, 100, int(ss.threshold), key="man_cons_t")
    with c2:
        man_op = st.radio("Critério usado", [">=", "<"],
                          index=0 if ss.cons_operator == ">=" else 1,
                          format_func=lambda op: "Conservados (≥)" if op == ">=" else "Únicos (<)",
                          horizontal=True, key="man_cons_op")
    if st.button("Registrar", key="cons_declare_btn"):
        conservancy_parser.write_metadata(csv_dir, man_t, man_op, source="manual")
        ss.threshold = int(man_t)
        ss.cons_operator = man_op
        st.rerun()


def step_conservancy():
    st.header("4 · Epitope Conservancy Analysis")
    if ss.poa1_result is None:
        st.warning("Execute o POA1 (etapa 3) primeiro.")
        return
    st.caption("O IEDB não oferece API nem versão local desta ferramenta. O POA reimplementa o "
               "cálculo localmente (Bui et al., 2007). Você também pode enviar os CSVs manualmente.")

    c1, c2 = st.columns([1, 2])
    with c1:
        ss.threshold = st.number_input("Limiar de identidade de sequência (%)", 1, 100, ss.threshold)
    with c2:
        ss.cons_operator = st.radio(
            "Critério", [">=", "<"],
            index=0 if ss.cons_operator == ">=" else 1,
            format_func=lambda op: "Conservados (≥ limiar)" if op == ">=" else "Únicos (< limiar)",
            horizontal=True)
    st.caption("O critério escolhido aqui decide quais proteínas contam como *match* e fica gravado "
               "junto dos CSVs. Ele precisa ser o mesmo objetivo da etapa 5 — o POA2 recusa a "
               "execução se divergirem, em vez de apenas rotular a coluna com o valor errado.")
    mode = st.radio("Modo", ["Reimplementação local (recomendado)", "Upload manual dos CSVs"])

    # POA1 wrote one epitope FASTA per species; each must be compared against its OWN protein set.
    epitope_files = sorted(ctx.conservancy_epitopes_dir.glob("*_epitopes.fasta"))
    species = [f.name[: -len("_epitopes.fasta")].upper() for f in epitope_files]
    if not species:
        st.warning("O POA1 não gerou nenhum arquivo de epítopos por espécie. Verifique a etapa 3.")
        return
    st.info(f"{len(species)} espécie(s) a analisar: " + ", ".join(species))

    if mode.startswith("Reimplementação"):
        st.write("Conjunto de proteínas para comparação. O padrão é o FASTA da etapa 1, "
                 "**filtrado por espécie**: os epítopos de cada espécie são comparados apenas com "
                 "as proteínas dela. Envie abaixo um conjunto próprio (ex.: o *world set* de "
                 "diversidade) para qualquer espécie que precise de um.")
        ref_counts = backend.count_proteins_by_specie(ss.proteins_path) if ss.proteins_path else {}
        for sp in species:
            up = st.file_uploader(f"FASTA de comparação para {sp} (opcional)",
                                  type=["fasta", "fa", "faa", "txt"], key=f"comp_{sp}")
            if up is not None:
                ss.comparison_sets[sp] = backend.save_upload(
                    ctx.inputs_dir, f"comparison_{sp}.fasta", up.getvalue())
            chosen = ss.comparison_sets.get(sp)
            if chosen:
                n_comp = len(backend.parse_fasta_headers(chosen))
                st.caption(f"↳ {sp}: conjunto próprio (`{Path(chosen).name}`), {n_comp} proteína(s).")
            else:
                n_comp = ref_counts.get(sp, 0)
                if n_comp:
                    st.caption(f"↳ {sp}: {n_comp} proteína(s) do -f.")
                else:
                    st.warning(f"↳ {sp}: nenhuma proteína dessa espécie no -f — será comparada com "
                               "o conjunto completo (identidades misturam espécies). Envie um "
                               "conjunto próprio ou corrija o FASTA da etapa 1.")
            # One protein cannot produce a conservancy: the epitope was predicted on it, so it
            # matches at 100% and the result is 100% for any threshold.
            if n_comp == 1:
                st.warning(f"↳ {sp}: com **uma única proteína** de comparação o limiar não "
                           "discrimina nada — todo epítopo sai com 100% de conservância, qualquer "
                           "que seja o valor escolhido acima. Envie um conjunto de diversidade "
                           "(ex.: o *world set* da espécie) para que o limiar tenha efeito.")

        if st.button("▶️ Calcular conservância (local)", type="primary"):
            with st.status("Calculando conservância…", expanded=True) as status:
                try:
                    written = conservancy_client.run_conservancy_for_dir(
                        str(ctx.conservancy_epitopes_dir), ss.proteins_path, float(ss.threshold),
                        str(ctx.conservancy_csv_dir),
                        specie_proteins={k: v for k, v in ss.comparison_sets.items() if k in species},
                        operator=ss.cons_operator)
                    ss.conservancy_ready = len(written) > 0
                    status.update(label=f"{len(written)} CSV(s) gerado(s).", state="complete")
                except Exception as exc:
                    status.update(label="Falha no cálculo de conservância.", state="error")
                    st.exception(exc)
    else:
        st.caption("O POA2 confere com que limiar os CSVs foram gerados antes de filtrar. Quando o "
                   "cabeçalho do arquivo traz essa informação (é o caso dos downloads do IEDB), ela "
                   "é lida daí; senão, informe abaixo.")
        ups = st.file_uploader("CSVs da Conservancy Analysis (IEDB)", type=["csv"],
                               accept_multiple_files=True, key="cons_csvs")
        if ups:
            # Any metadata already here describes the CSVs being replaced, not these.
            stale = ctx.conservancy_csv_dir / conservancy_parser.METADATA_FILENAME
            if stale.exists():
                stale.unlink()
            for up in ups:
                backend.save_upload(ctx.conservancy_csv_dir, up.name, up.getvalue())
            ss.conservancy_ready = True
            st.success(f"{len(ups)} CSV(s) carregado(s).")

        if list(ctx.conservancy_csv_dir.glob("*.csv")):
            _declare_uploaded_threshold()

    csvs = sorted(ctx.conservancy_csv_dir.glob("*.csv"))
    if csvs:
        st.divider()
        st.subheader("Resultado da conservância")
        overview = backend.conservancy_overview(ctx.conservancy_csv_dir)
        for error in overview.attrs.get("errors", []):
            st.warning(f"Não foi possível resumir o CSV: {error}")
        if not overview.empty:
            st.dataframe(overview, use_container_width=True, hide_index=True)
            st.caption("**Com match** = epítopos com ao menos uma proteína casando no limiar com "
                       "que o arquivo foi gerado. É sobre essa coluna que o filtro `-m` do POA2 "
                       "age na etapa 5.")
            silent = overview[overview["Com match"] == 0]["Espécie"].tolist()
            if silent:
                st.warning("Espécie(s) sem nenhum epítopo com match: " + ", ".join(silent) +
                           ". Elas serão excluídas se o filtro de proteínas com match (-m) for maior que zero.")
            st.plotly_chart(
                px.bar(overview.melt(id_vars="Espécie", value_vars=["Epítopos", "Com match"],
                                     var_name="Medida", value_name="count"),
                       x="Espécie", y="count", color="Medida", barmode="group",
                       title="Epítopos por espécie e quantos têm match"),
                use_container_width=True)

        for c in csvs:
            with st.expander(f"📄 {c.name}"):
                try:
                    st.dataframe(pd.read_csv(c), use_container_width=True, height=300)
                except (OSError, ValueError, pd.errors.ParserError) as exc:
                    st.warning(f"Não foi possível ler {c.name}: {exc}")
                _download(c, f"⬇️ {c.name}", "text/csv")

        st.caption("Gravados em:")
        st.code(str(ctx.conservancy_csv_dir), language=None)


# =========================================================================== STEP 5
def _resume_from_results():
    """
    Continue an analysis from what is already in ``results/``.

    The session state lives in the Streamlit process, but the files do not — so a run started in
    one environment (Windows) can be finished in another (WSL, the only place pyTMHMM builds)
    without redoing steps 1-4.
    """
    found = backend.resumable_run(ctx)
    csvs, fastas = found["csvs"], found["references"]
    if not csvs:
        return
    st.divider()
    st.info(f"Encontrei {len(csvs)} CSV(s) de conservância já gravados nesta análise "
            f"(`{ctx.conservancy_csv_dir.name}/`): " + ", ".join(c.name for c in csvs))
    st.caption("Use isto para continuar uma análise iniciada em outra sessão — por exemplo ao "
               "trocar do Windows para a WSL só para rodar o POA2. Abra a mesma análise dos dois "
               "lados: os arquivos estão na pasta, não na sessão.")

    if not fastas:
        st.warning(f"Nenhum FASTA de proteínas em `{ctx.inputs_dir.name}/` — o POA2 precisa do `-f`.")
        return
    # Label each candidate with the species its headers carry, so a file that merged
    # non-conforming headers is visible as such instead of being picked by accident.
    wanted = sorted({c.name[: -len("_conservancy.csv")].upper() for c in csvs})
    labels, paths = [], []
    for f in fastas:
        species = backend.species_in_fasta(str(f))
        covers = [sp for sp in wanted if sp in species]
        mark = "✅" if len(covers) == len(wanted) else ("⚠️" if covers else "❌")
        labels.append(f"{mark} {f.name} — espécies: {', '.join(species) or 'nenhuma'}")
        paths.append(f)
    default = next((i for i, f in enumerate(fastas)
                    if set(wanted) <= set(backend.species_in_fasta(str(f)))), 0)
    pick = st.selectbox("FASTA de proteínas (-f) usado naquela análise", range(len(labels)),
                        index=default, format_func=lambda i: labels[i])
    st.caption(f"O POA2 precisa das proteínas de: {', '.join(wanted)}. Um `-f` com espécies a mais "
               "(ex.: cabeçalhos não convertidos) faz um mesmo epítopo casar com vários registros.")
    if st.button("Retomar com estes arquivos"):
        ss.proteins_path = str(paths[pick])
        ss.conservancy_ready = True
        # also rebuild the POA1 table from the epitope FASTAs, so the results (etapa 6) and the
        # 2D/3D views (etapa 7) work in the resumed session too
        if ss.poa1_result is None:
            ss.poa1_result = backend.poa1_result_from_disk(ctx)
        st.rerun()


def step_poa2():
    st.header("5 · POA2 — conservação e topologia de membrana")
    if not ss.conservancy_ready:
        st.warning("Prepare a Conservancy Analysis (etapa 4) primeiro.")
        _resume_from_results()
        return
    if not ss.proteins_path:
        st.warning("O POA2 precisa do FASTA de proteínas (-f) — volte à etapa 1.")
        return
    # pyTMHMM gates only the *run*, not the step: everything below except the TMHMM pass is read
    # from the conservancy CSVs, so on Windows (where pyTMHMM does not build) the parameters and
    # their effect can still be tuned here before switching to the WSL side to execute.
    tmhmm_ok = tmhmm_client.is_available()

    # The threshold POA2 applies comes from the data, not from the form: it used to only rename a
    # column, so a run could report "conservados ≥70%" while filtering CSVs computed at ≤100%.
    source = conservancy_parser.inspect_directory(str(ctx.conservancy_csv_dir))
    if source.origin == "unknown":
        st.warning("Não foi possível descobrir com que limiar os CSVs desta pasta foram gerados — "
                   "não há metadados e o cabeçalho não traz o valor. O POA2 vai usar o limiar da "
                   "etapa 4 **sem conseguir validá-lo**; confira antes de usar os resultados.")
        t_value, op_value = int(ss.threshold), ss.cons_operator
    else:
        t_value, op_value = int(source.threshold), source.operator
        origem = ("metadados gravados na etapa 4" if source.origin == "meta"
                  else "cabeçalho dos próprios CSVs")
        st.info(f"CSVs gerados com limiar "
                f"**{conservancy_parser.format_threshold(source.threshold)}%** e critério "
                f"**{op_value}** (lido dos {origem}). São esses os valores que o POA2 aplica.")
        if source.meta and source.meta.species:
            weak = sorted(sp for sp, info in source.meta.species.items()
                          if not info.get("discriminating", True))
            if weak:
                st.warning("Espécie(s) cujo conjunto de comparação tinha uma única proteína — o "
                           "limiar não discriminou nada e a conservância saiu 100%: "
                           + ", ".join(weak) + ". Refaça a etapa 4 com um conjunto de diversidade "
                           "se estes resultados forem usados.")

    c1, c2, c3 = st.columns(3)
    with c1:
        objective = st.radio("Objetivo", ["conserved", "unique"],
                             index=1 if conservancy_parser.operator_family(op_value) == "lt" else 0,
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
        st.metric("Limiar aplicado (t)", f"{op_value} {t_value}%")
        identity_filter = st.checkbox(
            "Usar o limiar também como filtro", value=False,
            help="Por padrão o limiar só define a coluna de % de matches, como no POA original. "
                 "Marcado, o POA2 também descarta epítopos cuja identidade não satisfaça o "
                 "critério: ≥ exige que a identidade mínima alcance o limiar (conservado em todo o "
                 "conjunto); < exige que a máxima fique abaixo dele (único).")

    requested_symbol = ">=" if objective == "conserved" else "<"
    mismatch = (source.origin != "unknown"
                and conservancy_parser.operator_family(requested_symbol) != source.family)
    if mismatch:
        st.error(f"O objetivo escolhido (**{requested_symbol}**) não é o mesmo com que os CSVs "
                 f"foram gerados (**{source.operator}**). Volte à etapa 4 e recalcule com este "
                 "critério, ou escolha o objetivo correspondente. Filtrar com um critério e "
                 "rotular o resultado com outro é exatamente o que esta checagem impede.")

    params = dict(objective=objective, t=t_value, imin=imin, imax=imax, m=mmatch, rf=rf,
                  idf=identity_filter)

    # Where the epitopes are lost, computed from the CSVs with the parameters currently on screen.
    # POA2 only ever reported the survivors, so an empty result could not be told from a bound set
    # too tight — and this answers it before paying for a TMHMM run.
    funnel = backend.conservancy_funnel(ctx.conservancy_csv_dir, params)
    if funnel.empty:
        st.warning("Não foi possível calcular o efeito dos filtros. Confira os cabeçalhos e o "
                   "conteúdo dos CSVs na etapa 4; todos devem usar o mesmo limiar de identidade.")
    if not funnel.empty:
        st.subheader("Efeito dos filtros de conservância")
        kept, total = int(funnel["Epítopos"].iloc[-1]), int(funnel["Epítopos"].iloc[0])
        c1, c2 = st.columns([2, 3])
        with c1:
            st.dataframe(funnel, use_container_width=True, hide_index=True)
            st.metric("Epítopos que seguem para a topologia", f"{kept} de {total}")
        with c2:
            st.plotly_chart(px.bar(funnel, x="Epítopos", y="Etapa", orientation="h",
                                   title="Epítopos restantes após cada filtro"),
                            use_container_width=True)
        st.caption("Calculado a partir dos CSVs da etapa 4 com os parâmetros acima — ainda sem "
                   "rodar o TMHMM. A etapa seguinte do POA2 é a topologia de membrana, que não "
                   "descarta epítopos: apenas acrescenta as colunas `Portion_*`.")
        if kept == 0:
            st.error("Nenhum epítopo sobrevive a estes filtros — o POA2 vai produzir uma planilha "
                     "vazia. Afrouxe o filtro indicado acima como responsável pela queda.")

    if not tmhmm_ok:
        st.error("**pyTMHMM não está instalado neste ambiente** — sem ele o POA2 não roda. "
                 "Os parâmetros e o efeito deles acima continuam valendo: ajuste-os aqui e "
                 "execute do outro lado.")
        if sys.platform == "win32":
            st.warning(
                "No Windows, `pip install pyTMHMM` **não funciona**: não existe *wheel* para "
                "Windows, o build precisa do Microsoft C++ Build Tools e o pacote (1.3.6) não "
                "compila contra numpy 2.x (usa `np.int_t`). Rode a interface pela **WSL**, onde o "
                "ambiente já está pronto — veja `TESTING_WSL.md`:")
            st.code("wsl\nconda activate poa\ncd /mnt/c/Users/<você>/…/POA_Project\n"
                    "streamlit run poa/gui/app.py", language="bash")
            st.caption("Abra `http://localhost:8501` no navegador do Windows e **reabra esta mesma "
                       "análise** — a pasta `results/` é a mesma nos dois lados, então tudo o que "
                       "você já gerou continua valendo; basta ir direto para a etapa 5.")
        else:
            st.info("Instale no ambiente atual (precisa de compilador C e `numpy<2`):")
            st.code('pip install "numpy<2"\npip install --no-build-isolation pyTMHMM', language="bash")

    if st.button("▶️ Executar POA2", type="primary", disabled=mismatch or not tmhmm_ok):
        with st.status("Executando POA2 (conservância + TMHMM)…", expanded=True) as status:
            try:
                ss.poa2_result = backend.run_poa2(ctx, params, ss.proteins_path)
                status.update(label="POA2 concluído.", state="complete")
            except Exception as exc:
                status.update(label="POA2 falhou.", state="error")
                st.exception(exc)

    if ss.poa2_result is not None:
        st.divider()
        st.subheader("Saída do POA2")
        result = ss.poa2_result.results
        st.write(f"**{len(result)} epítopo(s) selecionado(s).** Critério aplicado: "
                 f"`{ss.poa2_result.type_symbol}`.")
        st.dataframe(result, use_container_width=True, height=300)
        st.caption("Os downloads e o resumo por topologia estão na etapa 6.")


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
