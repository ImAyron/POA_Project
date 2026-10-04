"""UI-agnostic orchestration for the GUI: working dirs, argument building, stage runners.

Kept free of Streamlit imports so it can be unit-tested and reused. The Streamlit app
(:mod:`poa.gui.app`) is a thin layer that calls into this module.
"""
from __future__ import annotations

import json
import os
import re
import unicodedata
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from ..core import pipeline
from ..core.models import STANDARD_COLUMNS
from ..services import realdata_import
from ..services.cache import Cache


#: repository root — ``poa/gui/backend.py`` -> ``poa/gui`` -> ``poa`` -> repo
REPO_ROOT = Path(__file__).resolve().parents[2]


def default_results_root() -> Path:
    """
    Where the GUI keeps its analyses: ``<repo>/results``.

    A fixed, findable folder inside the project instead of a per-session ``%TEMP%`` directory the
    user cannot locate and Windows may clear. Override with ``POA_RESULTS_DIR`` if the results
    should live elsewhere (e.g. a shared drive).

    Each analysis lives in its own ``DDMMAAAA-NOME_DO_TESTE`` subfolder (see :func:`create_run`),
    so re-running never overwrites an earlier one. A tree written before that convention existed
    sits directly under this root and is still usable — :func:`list_runs` reports it as legacy.
    """
    return Path(os.environ.get("POA_RESULTS_DIR") or (REPO_ROOT / "results"))


@dataclass
class WorkContext:
    """
    Working directory tree for one analysis.

    ``cache_root`` points the prediction cache outside the analysis folder so several analyses
    share it: cache keys are content hashes, so a second run over the same sequences reuses the
    IEDB/BepiPred results instead of re-submitting them. Left unset, the cache stays inside
    ``root`` (the pre-run-folder layout).
    """

    root: Path
    cache_root: Optional[Path] = None

    def __post_init__(self):
        self.root = Path(self.root)
        self.cache_root = Path(self.cache_root) if self.cache_root else None
        for d in (self.inputs_dir, self.mhcii_dir, self.poa1_out_dir,
                  self.conservancy_csv_dir, self.poa2_out_dir, self.cache_dir):
            d.mkdir(parents=True, exist_ok=True)

    @property
    def inputs_dir(self) -> Path:
        return self.root / "inputs"

    @property
    def mhcii_dir(self) -> Path:
        return self.root / "mhcii"

    @property
    def poa1_out_dir(self) -> Path:
        return self.root / "poa1_out"

    @property
    def conservancy_epitopes_dir(self) -> Path:
        # where POA1 writes the per-species epitope FASTAs
        return self.poa1_out_dir / "Conservancy Analysis"

    @property
    def conservancy_csv_dir(self) -> Path:
        return self.root / "conservancy_csv"

    @property
    def poa2_out_dir(self) -> Path:
        return self.root / "poa2_out"

    @property
    def cache_dir(self) -> Path:
        return self.cache_root or (self.root / "cache")

    def cache(self) -> Cache:
        return Cache(self.cache_dir)


# --------------------------------------------------------------------------- analysis folders
#: ``DDMMAAAA-NOME_DO_TESTE`` — the date first so the folders sort chronologically by name.
RUN_DIR_RE = re.compile(r"^(?P<date>\d{8})-(?P<label>.+)$")

#: Written into each analysis folder so the original (unslugged) name survives.
RUN_METADATA_FILENAME = "run.json"

#: Subfolders that mark a directory as a POA analysis tree.
_RUN_MARKERS = ("inputs", "poa1_out", "conservancy_csv", "poa2_out", "mhcii")


def slugify_run_name(name: str) -> str:
    """``'Teste DENV — world set' -> 'TESTE_DENV_WORLD_SET'`` (folder-safe, accent-free)."""
    text = unicodedata.normalize("NFKD", str(name))
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = re.sub(r"[^A-Za-z0-9]+", "_", text.upper()).strip("_")
    return text[:60] or "ANALISE"


def run_dir_name(name: str, when: Optional[date] = None) -> str:
    """The folder name for an analysis: ``DDMMAAAA-NOME_DO_TESTE``."""
    when = when or date.today()
    return f"{when:%d%m%Y}-{slugify_run_name(name)}"


@dataclass
class RunInfo:
    """One analysis folder under the results root, with enough state to choose between them."""

    path: Path
    name: str                   # folder name, or '' for the legacy tree
    label: str                  # the NOME_DO_TESTE part, or the name the user typed
    date: Optional[date]
    legacy: bool = False
    created_at: str = ""
    n_epitope_fastas: int = 0
    n_csvs: int = 0
    has_poa2: bool = False

    @property
    def progress(self) -> str:
        """Short human summary of how far this analysis got."""
        parts = []
        if self.n_epitope_fastas:
            parts.append(f"POA1: {self.n_epitope_fastas} espécie(s)")
        if self.n_csvs:
            parts.append(f"{self.n_csvs} CSV(s) de conservância")
        if self.has_poa2:
            parts.append("POA2 concluído")
        return " · ".join(parts) or "vazia"


def _is_run_tree(path: Path) -> bool:
    return any((path / marker).is_dir() for marker in _RUN_MARKERS)


def _describe_run(path: Path, name: str, label: str, when: Optional[date],
                  legacy: bool = False) -> RunInfo:
    created_at = ""
    meta_path = path / RUN_METADATA_FILENAME
    if meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            label = str(meta.get("label") or label)
            created_at = str(meta.get("created_at", ""))
        except (OSError, ValueError):
            pass
    epitopes = path / "poa1_out" / "Conservancy Analysis"
    return RunInfo(
        path=path,
        name=name,
        label=label,
        date=when,
        legacy=legacy,
        created_at=created_at,
        n_epitope_fastas=len(list(epitopes.glob("*_epitopes.fasta"))) if epitopes.is_dir() else 0,
        n_csvs=len(list((path / "conservancy_csv").glob("*.csv"))) if (path / "conservancy_csv").is_dir() else 0,
        has_poa2=any((path / "poa2_out").glob("*.xlsx")) if (path / "poa2_out").is_dir() else False,
    )


def list_runs(root: Path) -> List[RunInfo]:
    """
    Every analysis under ``root``, newest first.

    A tree written before the run-folder convention sits directly under ``root``; it is listed
    last, flagged ``legacy``, and keeps working where it is — nothing is moved.
    """
    root = Path(root)
    if not root.is_dir():
        return []

    runs: List[RunInfo] = []
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        match = RUN_DIR_RE.match(entry.name)
        if not match:
            continue
        try:
            when = datetime.strptime(match.group("date"), "%d%m%Y").date()
        except ValueError:
            continue
        runs.append(_describe_run(entry, entry.name, match.group("label"), when))

    runs.sort(key=lambda r: (r.date or date.min, r.name), reverse=True)
    if _is_run_tree(root):
        runs.append(_describe_run(root, "", "(análise anterior, fora de pasta datada)", None,
                                  legacy=True))
    return runs


def create_run(root: Path, name: str, when: Optional[date] = None) -> WorkContext:
    """
    Start a new analysis at ``root/DDMMAAAA-NOME_DO_TESTE`` and return its context.

    A name already used today reuses that folder rather than silently creating a near-duplicate —
    the point of the convention is that one test has one place.
    """
    root = Path(root)
    ctx = WorkContext(root / run_dir_name(name, when), cache_root=root / "cache")
    meta_path = ctx.root / RUN_METADATA_FILENAME
    if not meta_path.exists():
        meta_path.write_text(json.dumps({
            "label": str(name).strip(),
            "folder": ctx.root.name,
            "created_at": datetime.now().isoformat(timespec="seconds"),
        }, ensure_ascii=False, indent=2), encoding="utf-8")
    return ctx


def open_run(root: Path, path: Path) -> WorkContext:
    """Re-open an existing analysis (legacy tree included), sharing the root-level cache."""
    root, path = Path(root), Path(path)
    return WorkContext(path, cache_root=root / "cache")


# --------------------------------------------------------------------------- argument builders
def build_poa1_args(files: Dict[str, str], params: Dict, proteins_path: str, out_dir: str) -> SimpleNamespace:
    """
    Assemble a POA1 args namespace.

    Parameters:
        files: optional keys 'b2','b3','p','n','m','x' mapping to result file/dir paths.
        params: optional keys 'bmin','bmax','pmin','pmax','mhla','mic','xmin','xmax','e'.
        proteins_path: path to the protein FASTA (-f).
        out_dir: analysis output directory (-d).
    """
    return SimpleNamespace(
        b2=files.get("b2", ""),
        b3=files.get("b3", ""),
        bmin=int(params.get("bmin", 0)),
        bmax=int(params.get("bmax", 0)),
        p=files.get("p", ""),
        pmin=int(params.get("pmin", 0)),
        pmax=int(params.get("pmax", 0)),
        n=files.get("n", ""),
        m=files.get("m", ""),
        mhla=str(params.get("mhla", "DR")),
        mic=int(params.get("mic", 50)),
        x=files.get("x", ""),
        xmin=int(params.get("xmin", 0)),
        xmax=int(params.get("xmax", 0)),
        d=out_dir,
        f=proteins_path,
        e=str(params.get("e", "n")),
    )


def build_poa2_args(params: Dict, csv_dir: str, proteins_path: str, out_dir: str) -> SimpleNamespace:
    """
    Assemble a POA2 args namespace.

    Parameters:
        params: keys 'objective' ('conserved'|'unique'), 't', 'imin','imax','m','rf' and the
            threshold-handling switches 'idf' (apply the identity threshold as a filter, not just
            as the basis of the percent column) and 'strict' (refuse CSVs whose header declares a
            different threshold).
    """
    conserved = params.get("objective", "conserved") == "conserved"
    return SimpleNamespace(
        g="True" if conserved else None,
        l=None if conserved else "True",
        t=int(params["t"]),
        d=csv_dir,
        r=out_dir,
        imin=int(params.get("imin", 60)),
        imax=int(params.get("imax", 100)),
        m=int(params.get("m", 60)),
        f=proteins_path,
        rf=(None if params.get("rf", None) in (None, "") else int(params["rf"])),
        idf=bool(params.get("idf", False)),
        strict=bool(params.get("strict", False)),
    )


# --------------------------------------------------------------------------- stage runners
def run_poa1(ctx: WorkContext, files: Dict[str, str], params: Dict, proteins_path: str) -> pipeline.Poa1Result:
    args = build_poa1_args(files, params, proteins_path, str(ctx.poa1_out_dir))
    return pipeline.run_poa1(args)


def run_poa2(ctx: WorkContext, params: Dict, proteins_path: str) -> pipeline.Poa2Result:
    args = build_poa2_args(params, str(ctx.conservancy_csv_dir), proteins_path, str(ctx.poa2_out_dir))
    return pipeline.run_poa2(args)


# --------------------------------------------------------------------------- summaries for charts
def predictions_summary(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Return small aggregate tables (by method, by species) for charting the POA1 output."""
    out: Dict[str, pd.DataFrame] = {}
    if df is None or df.empty:
        return {"by_method": pd.DataFrame(columns=["Method", "count"]),
                "by_species": pd.DataFrame(columns=["Specie", "count"]),
                "length_hist": pd.DataFrame(columns=["length", "count"])}
    out["by_method"] = df.groupby("Method").size().reset_index(name="count")
    out["by_species"] = df.groupby("Specie").size().reset_index(name="count")
    lengths = df["Peptide Sequence"].astype(str).str.len()
    out["length_hist"] = lengths.value_counts().sort_index().reset_index()
    out["length_hist"].columns = ["length", "count"]
    return out


def topology_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate POA2 topology fractions (mean outer/TM/inner) for charting."""
    cols = ["Portion_Outside", "Portion_TM", "Portion_Inside"]
    if df is None or df.empty or not all(c in df.columns for c in cols):
        return pd.DataFrame(columns=["region", "mean_fraction"])
    means = {
        "Outside": pd.to_numeric(df["Portion_Outside"], errors="coerce").mean(),
        "Transmembrane": pd.to_numeric(df["Portion_TM"], errors="coerce").mean(),
        "Inside": pd.to_numeric(df["Portion_Inside"], errors="coerce").mean(),
    }
    return pd.DataFrame({"region": list(means.keys()), "mean_fraction": list(means.values())})


def save_upload(dest_dir: Path, filename: str, data: bytes) -> str:
    """Persist an uploaded file's bytes to ``dest_dir/filename`` and return the path."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / filename
    path.write_bytes(data)
    return str(path)


# --------------------------------------------------------------------------- multi-species merging
# POA1 takes ONE file per method (-b2/-b3/-p/-n/-x), but a multi-species run produces one file per
# species (one BepiPred submission per virus, one NetCTL page per virus, ...). These helpers merge
# the per-species files into the single input each method expects, so several species can be
# analysed in the same run instead of the last upload silently replacing the previous one.

def merge_bepipred2(paths: Sequence[str], out_path: str) -> str:
    """
    Merge several BepiPred-2.0 JSON exports into one, concatenating their ``antigens`` blocks.

    Each antigen key must already follow ``Protein_Specie_ID`` (the GUI adapter rewrites the
    generic ``Sequence`` key). Colliding keys get a ``#n`` suffix so nothing is dropped.
    """
    import json

    merged: Dict[str, dict] = {}
    info = {}
    for path in paths:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        info = data.get("info", info)
        for key, antigen in data.get("antigens", {}).items():
            new_key, n = key, 2
            while new_key in merged:
                new_key = f"{key}{n}"
                n += 1
            merged[new_key] = antigen
    if not merged:
        raise ValueError("No 'antigens' block found in any of the BepiPred-2.0 JSON files.")

    out = {"info": info, "antigens": merged} if info else {"antigens": merged}
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(out, fh)
    return out_path


def merge_fastas(paths: Sequence[str], out_path: str) -> str:
    """
    Concatenate FASTA files, keeping the first record for each duplicated header.

    Used both for the multi-species ``-f`` reference set and for merging per-species
    BepiPred-3.0 / 'other predictor' FASTAs.
    """
    from Bio import SeqIO

    seen = set()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as out:
        for path in paths:
            for rec in SeqIO.parse(path, "fasta"):
                if rec.id in seen:
                    continue
                seen.add(rec.id)
                out.write(f">{rec.description}\n{str(rec.seq)}\n")
    return out_path


def merge_text(paths: Sequence[str], out_path: str) -> str:
    """
    Concatenate plain-text prediction outputs (PAP/IMED ``.txt``, NetCTL ``.html``).

    Both parsers read the file line by line and take the species from each record's own header
    (or ``Protein_identifier`` column), so a plain concatenation is a valid multi-species input.
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as out:
        for path in paths:
            text = Path(path).read_text(encoding="utf-8", errors="replace")
            out.write(text)
            if not text.endswith("\n"):
                out.write("\n")
    return out_path


#: how each POA1 method key is merged when several per-species files are supplied
_MERGERS = {
    "b2": (merge_bepipred2, "bepipred2_merged.json"),
    "b3": (merge_fastas, "bepipred3_merged.fasta"),
    "p": (merge_text, "pap_merged.txt"),
    "n": (merge_text, "netctl_merged.html"),
    "x": (merge_text, "others_merged.fasta"),
}


def consolidate_sources(sources: Dict[str, List[str]], dest_dir: Path) -> Dict[str, str]:
    """
    Turn ``{method_key: [per-species files]}`` into the ``{method_key: single file}`` POA1 expects.

    ``'m'`` (MHC-II) is passed through untouched — it is already a directory of
    ``Protein_Specie.html`` files, which the parser walks per species.
    """
    files: Dict[str, str] = {}
    for key, paths in sources.items():
        paths = [p for p in paths if p]
        if not paths:
            continue
        if key == "m" or len(paths) == 1:
            files[key] = paths[0]
            continue
        merger, name = _MERGERS[key]
        files[key] = merger(paths, str(Path(dest_dir) / name))
    return files


#: residues that are not one of the 20 standard amino acids ('X' is the one POA1 refuses)
AMBIGUOUS_RESIDUES = frozenset("XBZJUO*")

#: how each POA1 method key is presented in the step-2 preview, and how its file is read
_SOURCE_KINDS: Dict[str, Tuple[str, str]] = {
    "b3": ("BepiPred-3.0", "fasta"),
    "b2": ("BepiPred-2.0", "json"),
    "p": ("PAP/IMED (antigenicidade)", "text"),
    "n": ("NetCTL 1.2", "text"),
    "m": ("MHC-II Binding (IEDB)", "dir"),
    "x": ("Outros preditores", "fasta"),
}


@dataclass
class SourcePreview:
    """
    What one prediction input holds, read from the file itself — no parser involved.

    Step 2 only collects raw tool output; POA1 (step 3) is what parses it, using filtering
    parameters that are chosen there. A preview that ran the parsers would therefore have to
    invent those parameters and could disagree with the real run. This reports what the *file*
    carries — records, headers, species — which is enough to catch the wrong file or a header
    convention the pipeline cannot read, without claiming to predict POA1's output.
    """
    label: str
    name: str
    kind: str
    size_bytes: int
    n_records: Optional[int] = None
    species: List[str] = field(default_factory=list)
    proteins: List[str] = field(default_factory=list)
    rows: List[Dict] = field(default_factory=list)
    head: str = ""
    notes: List[str] = field(default_factory=list)


def summarize_source(path: str, key: str) -> SourcePreview:
    """Describe one prediction file (see :class:`SourcePreview`); ``key`` is the POA1 method."""
    label, kind = _SOURCE_KINDS.get(key, (key, "text"))
    target = Path(path)
    preview = SourcePreview(label=label, name=target.name, kind=kind,
                            size_bytes=target.stat().st_size if target.exists() else 0)
    if not target.exists():
        preview.notes.append("Arquivo não encontrado.")
        return preview

    if kind == "dir" or target.is_dir():
        files = sorted(p for p in target.glob("*.htm*"))
        preview.kind = "dir"
        preview.n_records = len(files)
        preview.size_bytes = sum(p.stat().st_size for p in files)
        preview.rows = [{"Arquivo": p.name, "Bytes": p.stat().st_size} for p in files]
        if not files:
            preview.notes.append("Nenhum .html nesta pasta.")
        return preview

    if kind == "fasta":
        records = parse_fasta_headers(str(target))
        preview.n_records = len(records)
        preview.species = sorted({r["specie"] for r in records if r["specie"]})
        preview.proteins = sorted({r["protein"] for r in records if r["protein"]})
        preview.rows = [{"Cabeçalho": r["header"], "Proteína": r["protein"],
                         "Espécie": r["specie"], "Resíduos": r["length"]} for r in records]
        if records and not all(r["ok"] for r in records):
            preview.notes.append(
                "Há cabeçalhos sem `_`: nenhuma espécie sai deles. O POA1 agrupa esses epítopos "
                "sem espécie e a etapa 4 os compara com o conjunto errado.")
        if not records:
            preview.notes.append("Nenhum registro FASTA neste arquivo.")
        return preview

    if kind == "json":
        try:
            keys = realdata_import.antigen_keys(str(target))
            bad = nonconforming_antigen_keys(str(target))
        except Exception as exc:  # malformed JSON must not take the whole step down
            preview.notes.append(f"Não foi possível ler o JSON: {exc}")
            return preview
        preview.n_records = len(keys)
        preview.rows = [{"Antígeno": k, "Convenção": "❌" if k in bad else "✅"} for k in keys]
        preview.species = sorted({str(k).upper().split("_")[1] for k in keys
                                  if len(str(k).split("_")) > 1})
        if bad:
            preview.notes.append(
                "Antígeno(s) fora da convenção `Proteína_Espécie_ID`: " + ", ".join(map(str, bad))
                + ". Use o adaptador desta etapa para mapear cada um à sua espécie.")
        return preview

    text = target.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    preview.n_records = len(lines)
    preview.head = "\n".join(lines[:15])
    return preview


def summarize_sources(sources: Dict[str, List[str]]) -> List[SourcePreview]:
    """Describe every collected prediction input, in method order."""
    return [summarize_source(path, key)
            for key in sorted(sources) for path in sources[key]]


def conservancy_overview(csv_dir) -> pd.DataFrame:
    """
    One row per conservancy CSV: epitopes, how many matched, and the identity range.

    Step 4 used to list only file names, so the CSVs it had just produced could not be judged
    without downloading them. "Com match" counts epitopes with at least one protein match at the
    threshold the file was computed with — the quantity POA2's ``-m`` filter acts on.
    """
    from ..core.parsers import conservancy as cons

    def as_number(series):
        """Same '80.00%' -> 80.0 reading POA2 applies, so the summary cannot drift from it."""
        return cons.strip_percent(series)

    rows = []
    errors = []
    for path in cons.map_EpConservFiles(str(csv_dir)):
        try:
            df = pd.read_csv(path)
            percent_col = cons.find_percent_column(df.columns)
            minimum, maximum = as_number(df["Minimum identity"]), as_number(df["Maximum identity"])
        except (OSError, ValueError, KeyError, pd.errors.ParserError) as exc:
            errors.append(f"{Path(path).name}: {exc}")
            continue
        percent = as_number(df[percent_col].astype(str).str.split(" ").str[0])
        name = Path(path).name
        stem = (name[: -len("_conservancy.csv")] if name.endswith("_conservancy.csv")
                else Path(name).stem)
        rows.append({
            "Espécie": stem.upper(),
            "Arquivo": name,
            "Epítopos": len(df),
            "Com match": int((percent > 0).sum()),
            "Identidade mín. (%)": round(float(minimum.min()), 2) if len(df) else None,
            "Identidade máx. (%)": round(float(maximum.max()), 2) if len(df) else None,
        })
    result = pd.DataFrame(rows)
    result.attrs["errors"] = errors
    return result


#: PT-BR labels for the language-neutral stage keys of ``conservancy.stage_counts``
_FUNNEL_LABELS = {
    "input": lambda v: "Epítopos nos CSVs de conservância",
    "min_identity": lambda v: f"Identidade mínima ≥ {v}%",
    "max_identity": lambda v: f"Identidade máxima ≤ {v}%",
    "seq_match": lambda v: f"Proteínas com match ≥ {v}%",
    "identity_filter": lambda v: f"Identidade do epítopo {v[0]} {v[1]}%",
}


def conservancy_funnel(csv_dir, params: Dict) -> pd.DataFrame:
    """
    Where epitopes are lost along POA2's filters, as an ordered ``Etapa``/``Epítopos`` table.

    POA2 reports only the survivors, so an empty result cannot be told from a too-tight bound.
    ``params`` is the dict built by step 5 (``t``/``imin``/``imax``/``m``/``idf``/``objective``).
    """
    from ..core.parsers import conservancy as cons

    stages = cons.stage_counts(
        str(csv_dir),
        ID_threshold=params.get("t", 70),
        symbol=">=" if params.get("objective", "conserved") == "conserved" else "<",
        seq_match=params.get("m", 60),
        max_ID=params.get("imax", 100),
        min_ID=params.get("imin", 60),
        identity_filter=params.get("idf", False),
    )
    return pd.DataFrame([
        {"Etapa": _FUNNEL_LABELS[s["stage"]](s["value"]), "Epítopos": s["remaining"]}
        for s in stages
    ])


def parse_fasta_headers(fasta_path: str) -> List[Dict[str, str]]:
    """
    Show how each ``-f`` header is read under the ``Protein_Specie_ID`` convention.

    Returns one dict per record with ``header``/``protein``/``specie``/``id``/``ok``, plus
    ``length`` and ``ambiguous``. ``ok`` is False when the header has no ``_`` separator, so no
    species can be extracted at all — the case that makes POA1 report every epitope without a
    species. ``ambiguous`` lists the non-standard residues present: POA1's integrity check
    refuses an epitope containing ``X``, so showing it in step 1 explains a later refusal
    instead of letting it surface three steps on.
    """
    from Bio import SeqIO

    rows: List[Dict[str, str]] = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        parts = str(rec.id).upper().split("_")
        sequence = str(rec.seq).upper()
        rows.append({
            "header": str(rec.id),
            "protein": parts[0] if parts else "",
            "specie": parts[1] if len(parts) > 1 else "",
            "id": "_".join(parts[2:]) if len(parts) > 2 else "",
            "ok": len(parts) > 1,
            "length": len(sequence),
            "ambiguous": "".join(sorted(set(sequence) & AMBIGUOUS_RESIDUES)),
        })
    return rows


def normalize_fasta_headers(paths: Sequence[str], out_path: str,
                            swap: bool = False, default_id: str = "ref") -> str:
    """
    Rewrite FASTA headers to ``Protein_Specie_ID`` — sequences untouched.

    Parameters:
        swap: the source headers are ``Specie_Protein[_ID]`` (e.g. ``denv1_ns1``); swap the first
              two fields.
        default_id: appended when the header carries no third field.
    """
    from Bio import SeqIO

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as out:
        for path in paths:
            for rec in SeqIO.parse(path, "fasta"):
                parts = str(rec.id).upper().split("_")
                first = parts[0] if parts else str(rec.id).upper()
                second = parts[1] if len(parts) > 1 else ""
                rest = parts[2:] or [default_id]
                protein, specie = (second, first) if swap else (first, second)
                fields = [f for f in (protein, specie) if f] + rest
                out.write(">" + "_".join(fields) + "\n" + str(rec.seq) + "\n")
    return out_path


def species_in_fasta(fasta_path: str) -> List[str]:
    """Species tokens (field 2 of ``Protein_Specie_ID``) present in a protein FASTA."""
    species = []
    for row in parse_fasta_headers(fasta_path):
        if row["specie"] and row["specie"] not in species:
            species.append(row["specie"])
    return species


def count_proteins_by_specie(fasta_path: str) -> Dict[str, int]:
    """
    How many protein records each species contributes to a FASTA.

    Distinct from :func:`species_in_fasta`, which de-duplicates: the *size* of a species'
    comparison set is what decides whether a conservancy threshold can discriminate at all.
    """
    counts: Dict[str, int] = {}
    for row in parse_fasta_headers(fasta_path):
        if row["specie"]:
            counts[row["specie"]] = counts.get(row["specie"], 0) + 1
    return counts


def species_diagnostics(proteins_path: str, predictions: pd.DataFrame) -> Dict[str, List[str]]:
    """
    Compare the species found in the predictions with those in the ``-f`` reference set.

    Returns ``{'reference', 'predicted', 'missing_from_reference', 'unpredicted', 'unnamed'}``.
    ``missing_from_reference`` is the failure mode that silently produces meaningless conservancy
    results: epitopes of a species whose protein is not in ``-f``.
    """
    ref = species_in_fasta(proteins_path) if proteins_path else []
    if predictions is None or predictions.empty or "Specie" not in predictions:
        predicted, unnamed = [], []
    else:
        col = predictions["Specie"]
        unnamed = ["(sem espécie)"] if col.isna().any() else []
        predicted = sorted({str(v).upper() for v in col.dropna().unique()})
    return {
        "reference": ref,
        "predicted": predicted,
        "missing_from_reference": [sp for sp in predicted if sp not in ref],
        "unpredicted": [sp for sp in ref if sp not in predicted],
        "unnamed": unnamed,
    }


# --------------------------------------------------------------------------- real-data adapter
def import_bepipred2(ctx: WorkContext, filename: str, data: bytes,
                     specie: str, protein: str = "E",
                     mapping: Optional[Dict[str, Tuple[str, str, str]]] = None,
                     ) -> realdata_import.Prepared:
    """
    Adapt an original-workflow BepiPred-2.0 JSON upload to pipeline conventions.

    Persists the upload, then rewrites its antigen key(s) to ``Protein_Specie_ID`` and rebuilds the
    reference protein FASTA from the JSON's own ``AA`` array (see :mod:`poa.services.realdata_import`).
    Pass ``mapping`` (``{original_key: (protein, specie, id)}``) when one JSON holds antigens of
    different species — a single BepiPred run over several organisms.
    Returns the :class:`~poa.services.realdata_import.Prepared` record (paths + ``has_x`` flag);
    the caller wires ``b2_json`` into the POA1 files and may adopt ``reference_fasta`` as ``-f``.
    """
    raw_path = save_upload(ctx.inputs_dir, filename, data)
    return realdata_import.prepare_bepipred2(raw_path, specie, str(ctx.inputs_dir),
                                             protein=protein, mapping=mapping)


def predictions_from_epitope_fastas(ctx: "WorkContext") -> pd.DataFrame:
    """
    Rebuild the POA1 epitope table from the per-species FASTAs POA1 wrote.

    Their headers are ``Specie_Protein_Method_Init_Final`` and the sequence follows, which carries
    every standardized column except ``ID_Sequence``. Lets a resumed session show the results and
    the 2D/3D views without re-running POA1.
    """
    rows = []
    for fasta in sorted(ctx.conservancy_epitopes_dir.glob("*_epitopes.fasta")):
        header = None
        for line in fasta.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif line and header:
                fields = header.split("_")
                if len(fields) == 5:
                    specie, protein, method, init, final = fields
                    rows.append({"Method": method, "Specie": specie, "Protein": protein,
                                 "ID_Sequence": "-", "Initial Position": init,
                                 "Final Position": final, "Peptide Sequence": line})
                header = None
    return pd.DataFrame(rows, columns=list(STANDARD_COLUMNS))


def poa1_result_from_disk(ctx: "WorkContext") -> Optional[pipeline.Poa1Result]:
    """Reconstruct a :class:`~poa.core.pipeline.Poa1Result` from ``results/``; None if absent."""
    predictions = predictions_from_epitope_fastas(ctx)
    if predictions.empty:
        return None
    return pipeline.Poa1Result(
        predictions=predictions,
        report_path=str(ctx.poa1_out_dir / "Analysis_report.txt"),
        conservancy_dir=str(ctx.conservancy_epitopes_dir),
    )


def resumable_run(ctx: "WorkContext") -> Dict[str, List[Path]]:
    """
    What an earlier run left in the results tree, so a new session can pick it up.

    Session state lives in the Streamlit process but the files do not: a run started on Windows
    can be finished under WSL (the only environment where pyTMHMM builds) without repeating
    steps 1-4. Returns ``{'csvs': [...], 'references': [...]}`` — both empty when there is
    nothing to resume.
    """
    return {
        "csvs": sorted(ctx.conservancy_csv_dir.glob("*.csv")),
        "references": sorted(ctx.inputs_dir.glob("*.fasta")),
    }


def antigen_keys(path: str) -> List[str]:
    """Antigen keys of a saved BepiPred-2.0 JSON, for building the per-antigen mapping form."""
    return realdata_import.antigen_keys(path)


def nonconforming_antigen_keys(path: str) -> List[str]:
    """
    Antigen keys the BepiPred-2.0 parser cannot resolve into ``Protein_Specie_ID``.

    A key like ``Sequence`` or ``denv1`` yields a NaN species and the epitopes end up in a single
    ``nan_epitopes.fasta`` — checking up front lets the GUI send the file to the adapter instead
    of accepting it and failing two steps later.
    """
    import re

    from ..core.parsers.bepipred import BP2_ANTIGEN_KEY

    return [k for k in realdata_import.antigen_keys(path) if not re.match(BP2_ANTIGEN_KEY, str(k))]
