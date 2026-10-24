"""UI-agnostic orchestration for the GUI: working dirs, argument building, stage runners.

Kept free of Streamlit imports so it can be unit-tested and reused. The Streamlit app
(:mod:`poa.gui.app`) is a thin layer that calls into this module.
"""
from __future__ import annotations

import json
import hashlib
import os
import re
import shutil
import time
import unicodedata
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from ..core import pipeline
from ..core.models import STANDARD_COLUMNS
from ..core.artifacts import safe_child, file_hash
from ..logging_conf import TAG_WARN, get_logger
from ..services import realdata_import
from ..services.cache import Cache

logger = get_logger("gui")


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

    A name already used today gets a numeric suffix; creating never reopens an old run.
    """
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    base = run_dir_name(name, when)
    number = 1
    while True:
        path = root / (base if number == 1 else f"{base}_{number}")
        try:
            path.mkdir()
            break
        except FileExistsError:
            number += 1
    ctx = WorkContext(path, cache_root=root / "cache")
    # The loop above only breaks on a mkdir() that succeeded, so this folder was empty a moment
    # ago and its metadata is always ours to write — an existing run.json is never adopted.
    (ctx.root / RUN_METADATA_FILENAME).write_text(json.dumps({
        "label": str(name).strip(),
        "folder": ctx.root.name,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    return ctx


def open_run(root: Path, path: Path) -> WorkContext:
    """Re-open an existing analysis (legacy tree included), sharing the root-level cache."""
    root, path = Path(root), Path(path)
    return WorkContext(path, cache_root=root / "cache")


def delete_run(root: Path, path: Path) -> Dict[str, object]:
    """
    Permanently delete one analysis folder and everything inside it.

    This removes predictions, conservancy CSVs, POA2 spreadsheets and the run's log. There is no
    undo and nothing is moved to a recycle bin, so the caller must have asked first — the guards
    here only stop the mistakes a confirmation dialog cannot catch:

    * the path has to be **inside** ``root``, so a bad argument cannot walk out of the results
      tree and delete something else;
    * it must look like an analysis (``run.json`` or the usual subfolders), never an arbitrary
      directory someone happens to have under ``results/``;
    * ``root`` itself and the shared ``cache/`` are refused outright — the cache belongs to every
      analysis, not to this one.

    Args:
        root: the results root.
        path: the analysis folder to delete.

    Returns:
        ``{"folder", "files", "bytes"}`` describing what was removed.

    Raises:
        ValueError: the path is outside ``root``, is ``root``/``cache``, or is not an analysis.
        FileNotFoundError: there is nothing at ``path``.
    """
    root, path = Path(root).resolve(), Path(path).resolve()

    if not path.exists():
        raise FileNotFoundError(f"Não existe nada em {path}.")
    if not path.is_dir():
        raise ValueError(f"{path} não é uma pasta.")
    if path == root:
        raise ValueError("A pasta raiz de resultados não pode ser apagada — só uma análise "
                         "de cada vez.")
    if root not in path.parents:
        raise ValueError(f"{path} está fora de {root}; recusado por segurança.")
    if path.name == "cache":
        raise ValueError("A pasta 'cache' é compartilhada por todas as análises e não pertence "
                         "a nenhuma delas.")
    if not (path / RUN_METADATA_FILENAME).exists() and not _is_run_tree(path):
        raise ValueError(f"{path.name} não parece uma análise (sem {RUN_METADATA_FILENAME} nem "
                         f"as subpastas do pipeline). Apague-a pelo explorador de arquivos se "
                         f"for mesmo isso que você quer.")

    files = [p for p in path.rglob("*") if p.is_file()]
    total = sum(p.stat().st_size for p in files)
    summary = {"folder": path.name, "files": len(files), "bytes": total}

    shutil.rmtree(path)
    logger.info("analysis deleted: %s (%d file(s), %s)", path.name, len(files),
                human_size(total))
    return summary


def human_size(nbytes: float) -> str:
    """Byte count as a short human-readable string."""
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024 or unit == "GB":
            return f"{nbytes:.0f} {unit}" if unit == "B" else f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} GB"


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


#: ``str(path) -> (size, mtime_ns, digest)``. Streamlit re-runs the whole script on every widget
#: interaction and each run fingerprints the same inputs up to four times, so without this a
#: checkbox click costs a full SHA-256 pass over every prediction file — tens of MB each.
_DIGEST_CACHE: Dict[str, Tuple[int, int, str]] = {}

#: A digest is only remembered once the file has held still this long. Some mounts report mtime
#: at a coarse resolution, and a rewrite landing in the same tick as the hash would keep both the
#: size and the mtime — invisible to the check below. Freshly touched files are read again.
_DIGEST_SETTLE_NS = 2_000_000_000


def cached_file_hash(path) -> str:
    """:func:`file_hash` memoized on (size, mtime): the same digest without re-reading the bytes."""
    path = Path(path)
    stat = path.stat()
    key = str(path)
    remembered = _DIGEST_CACHE.get(key)
    if remembered is not None and remembered[:2] == (stat.st_size, stat.st_mtime_ns):
        return remembered[2]
    digest = file_hash(path)
    if time.time_ns() - stat.st_mtime_ns > _DIGEST_SETTLE_NS:
        _DIGEST_CACHE[key] = (stat.st_size, stat.st_mtime_ns, digest)
    return digest


def input_fingerprint(paths, params=None) -> str:
    """Content-based identity for a GUI stage, including missing/deleted inputs."""
    inputs = []
    for value in sorted(str(p) for p in paths if p):
        path = Path(value)
        if path.is_dir():
            files = sorted(p for p in path.rglob("*") if p.is_file()
                           and not any(part.startswith(".poa_") for part in p.relative_to(path).parts))
            inputs.append((value, [(str(p.relative_to(path)), cached_file_hash(p)) for p in files]))
        else:
            inputs.append((value, cached_file_hash(path) if path.is_file() else None))
    return hashlib.sha256(json.dumps([inputs, params], sort_keys=True, default=str).encode()).hexdigest()


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


def save_upload(dest_dir: Path, filename: str, data: bytes, *,
                ctx: Optional["WorkContext"] = None, role: str = "", method: str = "",
                specie: str = "", original_name: str = "") -> str:
    """
    Persist an uploaded file's bytes to ``dest_dir/filename`` and return the path.

    The name is validated before anything is written: an upload is named by whoever produced the
    file, and a name carrying a separator or a reserved device would land outside the analysis.

    Pass ``ctx`` together with ``role`` to record *where in the interface this file came from*
    (see :func:`record_input`). Without it the file is still written — the history is an
    addition, never a precondition for saving.
    """
    path = safe_child(dest_dir, filename)
    dest_dir.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    if ctx is not None and role:
        record_input(ctx, path, role=role, method=method, specie=specie,
                     original_name=original_name or filename)
    return str(path)


# --------------------------------------------------------------------------- input history
#: Written in the analysis folder: which step each file came from, so reopening the analysis can
#: put every one of them back where it was used instead of guessing from its name.
INPUT_HISTORY_FILENAME = "input_history.json"


def record_input(ctx: "WorkContext", path, *, role: str, method: str = "", specie: str = "",
                 original_name: str = "") -> None:
    """
    Record that ``path`` entered the analysis as ``role`` (optionally for ``method``/``specie``).

    People upload files by hand at several steps, and the file name alone does not say which step
    asked for it: ``others_0_x.fasta`` and a reference FASTA are both FASTAs in ``inputs/``, and a
    file the user named themselves carries no convention at all. Recording the provenance when the
    file arrives is the only account that stays true, so reopening an analysis can return each
    file to the exact step it was used in.

    Paths are stored **relative to the analysis folder**: the same folder is opened from Windows
    and from WSL, where its absolute path is different, and an absolute path recorded on one side
    would be meaningless on the other.

    Recording the same file as the same thing again writes nothing. Some of the save sites sit in
    code Streamlit reruns on every interaction — the merged ``-f`` is rebuilt whenever step 1
    renders — and rewriting the file there would churn ``recorded_at`` on every click of an
    analysis where nothing arrived.

    Failure to write the history never interrupts the run — the file itself is already saved, and
    losing a line of provenance is worth less than losing the upload.
    """
    entry = {
        "path": _relative_to_run(ctx, path),
        "role": role,
        "method": method,
        "specie": specie,
        "original_name": original_name or Path(path).name,
    }
    history = read_input_history(ctx)
    previous = next((item for item in history if item.get("path") == entry["path"]), None)
    if previous is not None and all(previous.get(k) == v for k, v in entry.items()):
        return
    entry["recorded_at"] = datetime.now().isoformat(timespec="seconds")
    history = [item for item in history if item.get("path") != entry["path"]]
    history.append(entry)
    target = ctx.root / INPUT_HISTORY_FILENAME
    try:
        target.write_text(json.dumps(history, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError as exc:
        logger.warning("%s gui/history | could not record %s in %s: %s",
                       TAG_WARN, entry["path"], INPUT_HISTORY_FILENAME, exc)


def read_input_history(ctx: "WorkContext") -> List[Dict[str, str]]:
    """The recorded provenance of this analysis's files; ``[]`` when there is none or it is unreadable."""
    target = ctx.root / INPUT_HISTORY_FILENAME
    if not target.is_file():
        return []
    try:
        data = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        logger.warning("%s gui/history | ignoring an unreadable %s: %s",
                       TAG_WARN, INPUT_HISTORY_FILENAME, exc)
        return []
    return [item for item in data if isinstance(item, dict) and item.get("path")]


def _relative_to_run(ctx: "WorkContext", path) -> str:
    """``<run>/inputs/x.fasta`` -> ``inputs/x.fasta``; an outside path is kept as given."""
    try:
        return Path(path).resolve().relative_to(Path(ctx.root).resolve()).as_posix()
    except ValueError:
        return str(path)


# --------------------------------------------------------------------------- step state
#: Written in the analysis folder beside the input history: what each step *held* when it was last
#: rendered — which files it had selected and which numbers were in its form.
STEP_STATE_FILENAME = "step_state.json"


def record_step_state(ctx: "WorkContext", **fields) -> None:
    """
    Record what a step currently holds, replacing the fields given and leaving the rest.

    :func:`record_input` answers a different question: it records that a file *arrived*, and that
    is a fact about the file which never stops being true. Several steps, though, hold a
    *selection* rather than a growing pile — step 2's uploaders replace the list for their method,
    BepiPred-2.0 and 3.0 evict each other because they are the same B-cell method, and "Limpar
    espécies importadas" drops the adapted imports with no file arriving at all. None of those
    choices leave a trace in a log of arrivals, so an analysis rebuilt from arrivals alone hands
    POA1 every file that ever passed through step 2 instead of the ones the person settled on.

    The form parameters are here for the same reason: a reopened analysis that offers back the
    MHC-II alleles' default while showing the HTML files produced with other alleles invites
    re-running the step with inputs that do not match what is on screen.

    Paths are stored relative to the analysis folder, as in :func:`record_input`, because the same
    folder is opened from Windows and from WSL. Writing fails silently for the same reason too:
    the files are already saved, and a lost preference is worth less than a broken step.
    """
    state = read_step_state(ctx)
    merged = dict(state)
    for key, value in fields.items():
        merged[key] = _relativise(ctx, value)
    if merged == state:
        return
    try:
        (ctx.root / STEP_STATE_FILENAME).write_text(
            json.dumps(merged, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError as exc:
        logger.warning("%s gui/history | could not record the step state in %s: %s",
                       TAG_WARN, STEP_STATE_FILENAME, exc)


def _relativise(ctx: "WorkContext", value):
    """
    Rewrite any path inside a recorded value as relative to the analysis folder.

    No ``isabs`` guard: ``POA_RESULTS_DIR`` is a documented user override and may be relative, and
    then the paths handed to this function are relative to the *working directory*, not to the
    analysis folder. Skipping them left the step state holding a CWD-relative path that
    ``_apply_step_state`` resolved as ``ctx.root / p``, duplicating the prefix — the selection was
    then dropped as missing and the step reopened empty with its file still on disk.

    Sending every string through :func:`_relative_to_run` is safe because it resolves both sides
    and returns the string untouched when it does not land inside the analysis folder, which is
    what happens to the form values recorded alongside the paths.
    """
    if isinstance(value, str):
        return _relative_to_run(ctx, value)
    if isinstance(value, dict):
        return {k: _relativise(ctx, v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_relativise(ctx, v) for v in value]
    return value


def read_step_state(ctx: "WorkContext") -> Dict[str, object]:
    """What each step held when it was last rendered; ``{}`` when nothing was recorded."""
    target = ctx.root / STEP_STATE_FILENAME
    if not target.is_file():
        return {}
    try:
        data = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        logger.warning("%s gui/history | ignoring an unreadable %s: %s",
                       TAG_WARN, STEP_STATE_FILENAME, exc)
        return {}
    return data if isinstance(data, dict) else {}


#: How a file in ``inputs/`` is classified when nothing recorded where it came from. Matched in
#: order against the name, first hit wins. Each entry is the *stem* a save site in ``app.py``
#: writes ahead of the name the user gave; the separator is added at match time, so a stem never
#: matches the middle of a word. Both separators occur: '_' where the save site prefixes a name
#: the person chose ('ref_0_minhas.fasta'), '.' where it writes the whole name itself
#: ('bepipred3.fasta', from the automated run). This is a fallback, not the contract: a file
#: recorded by :func:`record_input` is read back from that record instead. It exists because files
#: already on disk have no history to read.
_LEGACY_INPUT_PATTERNS: Tuple[Tuple[str, str, str], ...] = (
    # (stem written by app.py, role, method)
    ("proteins_all", "merged_reference", ""),
    ("proteins_normalized", "merged_reference", ""),
    ("ref", "reference", ""),
    ("comparison", "comparison_set", ""),
    ("structure", "structure", ""),
    ("bepipred3", "prediction", "b3"),
    ("pap", "prediction", "p"),
    ("netctl", "prediction", "n"),
    ("others", "prediction", "x"),
    ("raw", "bepipred2_raw", ""),
)

#: The two files the adapter (:func:`poa.services.realdata_import.prepare_bepipred2`) writes. It
#: names both ``<stem><suffix>`` where the stem is ``'-'.join(species)`` — the species tokens
#: somebody typed into the step-2 mapping form.
_ADAPTER_SUFFIXES: Tuple[Tuple[str, str, str], ...] = (
    ("_b2.json", "bepipred2_adapted", "b2"),
    ("_ref.fasta", "bepipred2_reference", ""),
)


def _classify_legacy_input(name: str) -> Tuple[str, str]:
    """``(role, method)`` for a file in ``inputs/`` whose arrival was never recorded."""
    lowered = name.lower()
    # Both accounts can claim the same name, and which one is right is decided by what sits in
    # front of the suffix. A save site writes '<stem>_<a name the person chose>', and that name
    # may itself end in '_ref' — 'ref_0_E_DENV1_ref.fasta' is the step-1 upload the help text
    # suggests, not a rebuilt reference. The adapter's stem is a species token, so it carries no
    # '_' of its own: that absence is what separates 'DENV1_ref.fasta' from 'ref_0_E_DENV1_ref.
    # fasta', and it is also why 'PAP_b2.json' is the adapter's and not an antigenicity upload.
    for suffix, role, method in _ADAPTER_SUFFIXES:
        if lowered.endswith(suffix) and "_" not in lowered[: -len(suffix)]:
            return role, method
    for stem, role, method in _LEGACY_INPUT_PATTERNS:
        if lowered.startswith((f"{stem}_", f"{stem}.")):
            return role, method
    # POA0 writes '<PROT>_<SPECIE>.fasta' with no prefix of its own, and a hand-named reference
    # looks the same; both are references, which is the only role an unprefixed FASTA has here.
    if lowered.endswith((".fasta", ".fa", ".faa")):
        return "reference", ""
    return "", ""


def _is_assembled_input(name: str) -> bool:
    """
    True for a file :func:`consolidate_sources` assembles out of the inputs of one method.

    Nobody supplies these: they are the single file per method POA1 takes, rebuilt from the
    per-species parts on every run. Their names carry the same stem as the parts
    (``bepipred3_merged.fasta`` next to ``bepipred3_0_a.fasta``), so classifying one by name would
    offer POA1 the parts *and* their merge, counting every epitope in them twice. The names are
    read from the table that writes them, so the two cannot drift apart.
    """
    return name in {merged for _, merged in _MERGERS.values()}


def _usable_as_b2(path: Path) -> bool:
    """
    Whether a raw BepiPred-2.0 export is the ``-b2`` input itself, with no adapter in between.

    Step 2 asks this of every upload: keys already in ``Protein_Specie_ID`` are used as they are,
    generic ones (``Sequence``, ``denv1``) have to be mapped first and it is the adapter's output
    that the analysis then uses. The answer is a property of the file, not of the session, so a
    reopened analysis asks it again instead of trusting a recorded verdict — which is also what
    reaches the files that were on disk before any provenance was recorded.

    Handing an unmapped file back would be worse than handing nothing: POA1 resolves no species
    from it and every epitope lands under ``nan``.
    """
    try:
        return bool(antigen_keys(str(path))) and not nonconforming_antigen_keys(str(path))
    except Exception as exc:  # noqa: BLE001 - a malformed file must not break the reopening
        logger.warning("%s gui/history | %s: unreadable as BepiPred-2.0, so step 2 will ask for "
                       "it again: %s", TAG_WARN, path.name, exc)
        return False


@dataclass
class RestoredAnalysis:
    """Everything a reopened analysis could put back, and where it came from."""

    references: List[str] = field(default_factory=list)      # step 1, as uploaded/generated
    merged_reference: Optional[str] = None                   # step 1, the single -f
    prediction_files: Dict[str, List[str]] = field(default_factory=dict)   # step 2, by method
    prepared_b2: List[object] = field(default_factory=list)   # step 2, adapter records
    comparison_sets: Dict[str, str] = field(default_factory=dict)          # step 4, by species
    structure: Optional[str] = None                           # step 7
    conservancy_csvs: List[str] = field(default_factory=list)  # step 4
    poa1_result: Optional[object] = None                      # step 3
    poa2_result: Optional[object] = None                      # step 5, shown by step 6
    ref_swap: bool = False                                    # step 1, headers were Specie_Protein
    params: Dict[str, object] = field(default_factory=dict)   # the forms, by step
    source: str = "none"              # 'history' | 'inferred' | 'mixed' (both accounts) | 'none'

    @property
    def anything(self) -> bool:
        # Whatever summary() announces has to count here, because _restore_analysis_files gives up
        # when this is false: a folder whose only surviving artefact is a POA1 or POA2 result would
        # otherwise produce a summary nobody gets to see, and step 6 would keep saying POA2 never
        # ran — the exact case poa2_result_from_disk exists to cover, when the inputs are gone.
        return bool(self.references or self.merged_reference or self.prediction_files
                    or self.prepared_b2 or self.comparison_sets or self.structure
                    or self.conservancy_csvs
                    or self.poa1_result is not None or self.poa2_result is not None)

    def summary(self) -> List[str]:
        """Human lines describing what was restored, for the interface to show."""
        lines = []
        if self.merged_reference or self.references:
            lines.append(f"Etapa 1: {len(self.references) or 1} FASTA(s) de proteínas")
        for method, files in sorted(self.prediction_files.items()):
            lines.append(f"Etapa 2 ({method}): {len(files)} arquivo(s)")
        if self.prepared_b2:
            lines.append(f"Etapa 2 (b2): {len(self.prepared_b2)} JSON(s) adaptado(s)")
        if self.poa1_result is not None:
            lines.append(f"Etapa 3: {len(self.poa1_result.predictions)} epítopo(s)")
        if self.poa2_result is not None:
            lines.append(f"Etapa 5: {len(self.poa2_result.results)} epítopo(s) selecionado(s) "
                         f"({self.poa2_result.type_symbol})")
        if self.comparison_sets:
            lines.append("Etapa 4: conjunto de comparação de "
                         + ", ".join(sorted(self.comparison_sets)))
        if self.conservancy_csvs:
            lines.append(f"Etapa 4: {len(self.conservancy_csvs)} CSV(s) de conservância")
        if self.structure:
            lines.append(f"Etapa 7: estrutura {Path(self.structure).name}")
        return lines


def restore_analysis(ctx: "WorkContext") -> RestoredAnalysis:
    """
    Recover every input of an analysis folder, each tagged with the step it belongs to.

    Reopening an analysis used to restore only the conservancy CSVs and the reference, behind a
    button on step 5. Everything else — the prediction files people upload by hand at step 2, the
    per-species comparison sets of step 4, the structure of step 7 — stayed on disk while the
    interface showed empty fields, so continuing a run meant uploading it all again.

    Three accounts are read, in decreasing order of authority.

    :func:`record_step_state` is first where it speaks, because it is the only one that records a
    *choice*. Step 1 and step 2 replace what they hold rather than add to it, and the two BepiPred
    versions evict each other, so a step rebuilt from arrivals alone offers POA1 every file that
    ever passed through it — the JSON the person replaced, the BepiPred-2.0 import the automated
    3.0 run discarded — merged with the ones they settled on.

    :func:`record_input` is next: it says which step each file arrived at, which is what no file
    name can be trusted to say. Last comes the name (:func:`_classify_legacy_input`), best-effort,
    for files that were on disk before any of this was recorded; it only fills what the records
    leave empty, and ``source`` reports whenever it was used so the interface can ask for a check.
    """
    restored = RestoredAnalysis()
    history = read_input_history(ctx)
    state = read_step_state(ctx)
    recorded = {item["path"] for item in history}

    inferred = []
    for path in sorted(ctx.inputs_dir.glob("*")) if ctx.inputs_dir.is_dir() else []:
        if not path.is_file() or _relative_to_run(ctx, path) in recorded:
            continue
        if _is_assembled_input(path.name):
            continue
        role, method = _classify_legacy_input(path.name)
        if role:
            inferred.append((path, role, method, ""))

    recorded_entries = [(ctx.root / item["path"], item.get("role", ""), item.get("method", ""),
                         item.get("specie", "")) for item in history]
    restored.source = ("mixed" if inferred and history
                       else "inferred" if inferred or not history
                       else "history")

    adapted: Dict[str, Path] = {}

    def absorb(entries, *, authoritative: bool) -> None:
        """
        Place each file on its step. Guessed files fill only what the record left empty.

        Within the record the later entry wins — it is chronological, and the second structure
        somebody uploads is the one the analysis used. A guess never displaces a record, in a
        single slot or at the head of a list: ``merge_fastas`` keeps the first record of each
        header, so a stale legacy FASTA arriving first would hand POA1 the sequence the person
        reopened the analysis to replace.
        """
        for path, role, method, specie in entries:
            if not path.is_file():
                continue            # the file was deleted after being recorded
            text = str(path)
            if role == "reference":
                restored.references.append(text)
            elif role == "merged_reference":
                if authoritative or restored.merged_reference is None:
                    restored.merged_reference = text
            elif role == "prediction" and method == "m":
                # MHC-II is the one method POA1 takes as a *directory* of Protein_Specie.html
                # files. Each file is recorded on its own — that is the provenance — but the input
                # the step has to offer back is the folder, not the individual paths.
                pass
            elif role == "prediction" and method:
                restored.prediction_files.setdefault(method, []).append(text)
            elif role == "comparison_set":
                key = specie or path.stem.replace("comparison_", "").upper()
                if authoritative or key not in restored.comparison_sets:
                    restored.comparison_sets[key] = text
            elif role == "structure":
                if authoritative or restored.structure is None:
                    restored.structure = text
            elif role == "bepipred2_adapted":
                adapted.setdefault(path.name[: -len("_b2.json")], path)
            elif role == "bepipred2_raw" and _usable_as_b2(path):
                # A raw export whose keys already conform is the -b2 input itself: step 2 skips
                # the adapter and uses this very file, so this is the only record of it there is.
                restored.prediction_files.setdefault("b2", []).append(text)

    absorb(recorded_entries, authoritative=True)
    absorb(inferred, authoritative=False)

    # The adapter's two files are one record: pair them back up by their shared stem.
    for stem, b2_path in sorted(adapted.items()):
        prepared = realdata_import.prepared_from_disk(
            str(b2_path), str(ctx.inputs_dir / f"{stem}_ref.fasta"))
        if prepared is not None:
            restored.prepared_b2.append(prepared)

    # MHC-II is a directory of HTML files, not one upload, so it is read as such.
    if ctx.mhcii_dir.is_dir() and any(ctx.mhcii_dir.glob("*.htm*")):
        restored.prediction_files.setdefault("m", [str(ctx.mhcii_dir)])

    _apply_step_state(ctx, restored, state)

    if ctx.conservancy_csv_dir.is_dir():
        restored.conservancy_csvs = [str(p) for p in sorted(ctx.conservancy_csv_dir.glob("*.csv"))]

    restored.poa1_result = poa1_result_from_disk(ctx)
    restored.poa2_result = poa2_result_from_disk(ctx)
    return restored


def _apply_step_state(ctx: "WorkContext", restored: RestoredAnalysis,
                      state: Dict[str, object]) -> None:
    """
    Narrow what the files account for down to what the steps had actually selected.

    Only the fields that were recorded are applied: an analysis whose step state predates this —
    or that never rendered the step — keeps everything the files gave it, which is the best that
    can be said about it. A file the selection names but that is gone from disk is dropped, same
    as anywhere else in the restore.
    """
    def present(values) -> List[str]:
        paths = [str(ctx.root / p) for p in values if isinstance(p, str)]
        return [p for p in paths if Path(p).is_file()]

    if isinstance(state.get("references"), list):
        restored.references = present(state["references"])
    if isinstance(state.get("sources"), dict):
        selected: Dict[str, List[str]] = {}
        for method, values in state["sources"].items():
            if method == "m":
                continue
            files = present(values if isinstance(values, list) else [])
            if files:
                selected[method] = files
        # MHC-II is the folder itself, never a selection among files: the step hands POA1
        # ``mhcii/`` whenever it holds anything at all, so what is on disk already is the answer.
        if "m" in restored.prediction_files:
            selected["m"] = restored.prediction_files["m"]
        restored.prediction_files = selected
    if isinstance(state.get("prepared_b2"), list):
        restored.prepared_b2 = [pr for pr in restored.prepared_b2
                                if _relative_to_run(ctx, pr.b2_json) in state["prepared_b2"]]
    if isinstance(state.get("ref_swap"), bool):
        restored.ref_swap = state["ref_swap"]
    restored.params = {key: value for key, value in state.items()
                       if key.startswith("params_") and isinstance(value, dict)}


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
    Concatenate FASTAs, deduplicating identical records and rejecting conflicting IDs.

    Used both for the multi-species ``-f`` reference set and for merging per-species
    BepiPred-3.0 / 'other predictor' FASTAs.
    """
    from Bio import SeqIO

    seen = {}
    records = []
    for path in paths:
        for rec in SeqIO.parse(path, "fasta"):
            key = rec.id.upper()
            sequence = str(rec.seq)  # case encodes BepiPred-3 predictions; preserve it
            if key in seen:
                prior_sequence, prior_path = seen[key]
                if prior_sequence != sequence:
                    raise ValueError(f"Conflicting FASTA ID {rec.id!r}: {prior_path} and {path}")
                continue
            seen[key] = (sequence, path)
            records.append(rec)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as out:
        for rec in records:
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

    All three files this leaves in ``inputs/`` are recorded (:func:`record_input`). The adapter
    writes its two itself, so nothing else would record them, and the name scan can read them
    back — the pair is found by its shared stem either way. What the record buys is that they are
    *accounted for*: without it, reopening the analysis meets two files it has to guess at, and an
    analysis whose every file went through the interface is asked to check its steps for
    guesswork that never happened.
    """
    raw_path = save_upload(ctx.inputs_dir, filename, data, ctx=ctx, role="bepipred2_raw",
                           original_name=filename)
    prepared = realdata_import.prepare_bepipred2(raw_path, specie, str(ctx.inputs_dir),
                                                 protein=protein, mapping=mapping)
    # The species token is the stem both files share, and what pairs them back up. The reference
    # is recorded even though no branch of restore_analysis reads that role: it is there to say
    # the file is spoken for, so the scan does not meet an unexplained FASTA in inputs/.
    stem = "-".join(prepared.species)
    record_input(ctx, prepared.b2_json, role="bepipred2_adapted", method="b2", specie=stem)
    record_input(ctx, prepared.reference_fasta, role="bepipred2_reference", specie=stem)
    return prepared


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


def poa2_result_from_disk(ctx: "WorkContext") -> Optional[pipeline.Poa2Result]:
    """
    Reconstruct a :class:`~poa.core.pipeline.Poa2Result` from ``poa2_out/``; None if absent.

    The counterpart of :func:`poa1_result_from_disk`, and missing for as long as that one existed
    alone: reopening a finished analysis brought step 3 back but left step 6 saying the POA2 had
    not been run, with the spreadsheet it is asking for sitting in the folder. Running it again to
    read a file that is already there is the kind of work reopening exists to spare — and, with
    the filter form restored from :func:`read_step_state` only since this change, a re-run could
    quietly produce a *different* spreadsheet from the one lying beside it.

    :func:`poa.core.pipeline.run_poa2` writes ``POA2_analysis_<ge|lt><threshold>.xlsx``, so the
    name carries the criterion back; the optional ``-rf`` FASTA shares the stem.
    """
    if not ctx.poa2_out_dir.is_dir():
        return None
    sheets = sorted(ctx.poa2_out_dir.glob("POA2_analysis_*.xlsx"),
                    key=lambda p: p.stat().st_mtime)
    if not sheets:
        return None
    xlsx = sheets[-1]
    try:
        results = pd.read_excel(xlsx)
    except Exception as exc:  # noqa: BLE001 - a reopened analysis must not fail on a stale file
        logger.warning("%s gui/history | %s: could not be read back, so step 5 will ask for the "
                       "POA2 to be run again: %s", TAG_WARN, xlsx.name, exc)
        return None
    fastas = sorted(ctx.poa2_out_dir.glob(f"{xlsx.stem}_rf*.fasta"))
    return pipeline.Poa2Result(
        results=results,
        xlsx_path=str(xlsx),
        fasta_path=str(fastas[0]) if fastas else None,
        type_symbol=">=" if "_ge" in xlsx.stem else "<",
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


#: UniProt's own two accession patterns, as the registry defines them. Worth spelling out rather
#: than approximating: a first draft used ``[A-NR-Z]`` for the leading letter, which excludes O, P
#: and Q — and ``P29990`` (dengue NS1) is exactly the shape it rejected.
_UNIPROT_ACCESSION = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")


def accession_from_header(track: str, proteins_path: Optional[str]) -> str:
    """
    The UniProt accession of the reference behind ``Protein|Specie``, or ``''``.

    The convention puts it in the third field of the header — ``NS1_DENV1_P29990`` — which is what
    POA0 writes and what people working from UniProt tend to use as the ID. Anything else is not
    an accession and is better treated as absent than guessed at: PDBe answers by accession, and
    a wrong one would list the structures of another protein entirely, with nothing on screen to
    say so.
    """
    if not proteins_path:
        return ""
    protein, specie = track.split("|")
    for row in parse_fasta_headers(proteins_path):
        if row["protein"] == protein.upper() and row["specie"] == specie.upper():
            candidate = str(row["id"]).split("_")[0].strip().upper()
            if _UNIPROT_ACCESSION.fullmatch(candidate):
                return candidate
            return ""
    return ""


#: Written by POA0 next to the files it generates.
POA0_REPORT_FILENAME = "POA0_relatorio.json"


def poa0_report(ctx: "WorkContext") -> Optional[Dict]:
    """The POA0 report of this analysis, or ``None`` when the reference did not come from POA0."""
    path = ctx.inputs_dir / POA0_REPORT_FILENAME
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        logger.warning("%s gui/poa0 | ignoring an unreadable %s: %s",
                       TAG_WARN, POA0_REPORT_FILENAME, exc)
        return None


def structure_mapping(ctx: "WorkContext", pdb_path: Optional[str]) -> Optional[Dict[str, object]]:
    """
    The SIFTS correspondence POA0 recorded for ``pdb_path``, when it is the structure POA0 fetched.

    SIFTS states, residue by residue, which position of the UniProt sequence a structure's
    numbering corresponds to. It is exact, where an alignment is an estimate — and POA0 already
    computes and stores it (``estrutura.offset_regiao``), relative to the region it wrote, which
    is the coordinate system the epitope positions use.

    Until now nothing read it back: ``map_epitopes`` accepted ``sifts_offset`` and no caller ever
    passed one, so every structure fell to the alignment path while step 1 told the user the
    mapping was in use. This is the function that makes that claim true.

    Returns ``{"offset", "span", "pdb_id", "covers_region"}``, or ``None`` when this analysis has
    no POA0 report, when the report is about a different structure file, or when PDBe knew no
    mapping — all cases where the alignment remains the right answer.
    """
    if not pdb_path:
        return None
    report = poa0_report(ctx)
    structure = (report or {}).get("estrutura") or {}
    recorded = (report or {}).get("arquivos", {}).get("estrutura")
    if not structure or "offset_regiao" not in structure:
        return None
    # The report describes one structure. A different file in the viewer is not it, and applying
    # another structure's offset would place the epitopes confidently on the wrong residues.
    if not recorded or Path(recorded).name != Path(pdb_path).name:
        return None
    sifts = structure.get("sifts") or {}
    region = (report or {}).get("regiao") or {}
    length = int((report or {}).get("sequencia", {}).get("tamanho", 0))
    span = (1, length) if length else None
    return {
        "offset": int(structure["offset_regiao"]),
        "span": span,
        "pdb_id": structure.get("pdb_id", ""),
        "chain": sifts.get("chain", ""),
        "covers_region": bool(structure.get("cobre_regiao", True)),
        "region": (region.get("start"), region.get("end")),
    }


def reference_sequence(ctx: "WorkContext", protein: str, specie: str,
                       proteins_path: Optional[str] = None) -> Optional[Dict[str, str]]:
    """
    Find the sequence the predictions for one protein/species ran on.

    Epitope positions are relative to this sequence, so it is what a structure has to be
    compared against before its residues can be colored. Searches the analysis' own ``inputs/``
    for a record whose header reads ``Protein_Specie_ID``.

    Parameters:
        ctx: the analysis folder.
        protein: field 0 of the header.
        specie: field 1 of the header.

    Returns:
        ``{"header", "sequence", "path"}`` for the first match, or None when no input FASTA
        carries that protein/species — in which case the structure cannot be verified and the
        caller must say so rather than highlight anything.
    """
    from Bio import SeqIO

    manifest_path = ctx.poa1_out_dir / "poa1_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        path = safe_child(ctx.poa1_out_dir, manifest["reference"])
        if file_hash(path) != manifest["reference_sha256"]:
            raise ValueError("A referência salva do POA1 foi alterada; execute o POA1 novamente.")
    elif proteins_path:
        path = Path(proteins_path)
    else:
        return None  # legacy runs need an explicitly selected reference
    if not path.is_file():
        return None
    matches = []
    for rec in SeqIO.parse(str(path), "fasta"):
        parts = str(rec.id).upper().split("_")
        if len(parts) > 1 and parts[:2] == [protein.upper(), specie.upper()]:
            matches.append({"header": str(rec.id), "sequence": str(rec.seq).upper(), "path": str(path)})
    if len({r["sequence"] for r in matches}) > 1:
        raise ValueError("Há referências diferentes para esta proteína/espécie; selecione uma referência inequívoca antes de destacar epítopos.")
    return matches[0] if matches else None
