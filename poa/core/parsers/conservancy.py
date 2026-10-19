"""Parser and filter for IEDB Epitope Conservancy Analysis results (CSV files).

The filtering logic is the original ``conservancyAnalysis.py``. What changed is that the
sequence-identity threshold and its operator no longer merely *label* a column: they are read back
from the data — from the CSV header itself, the way IEDB spells it
(``Percent of protein sequence matches at identity >= 70%``), or from the ``conservancy_meta.json``
that :mod:`poa.services.conservancy_client` writes next to the CSVs — and checked against what
POA2 was asked to do (:func:`check_threshold`). Running POA2 with a ``-t``/``-g|-l`` that
contradicts the data being filtered is a silent scientific error, so it is now reported instead of
being renamed away.

The *local reimplementation* that generates such results lives in
:mod:`poa.services.conservancy_client`.
"""
from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import pandas as pd

from ...logging_conf import TAG_WARN, get_logger
from ..artifacts import file_hash, safe_child

logger = get_logger("conservancy")


#: Columns of an IEDB Epitope Conservancy Analysis CSV (the percent column carries the threshold).
IEDB_COLUMNS = [
    "Epitope #",
    "Epitope name",
    "Epitope sequence",
    "Epitope length",
    "Percent of protein sequence matches at identity <= 100%",
    "Minimum identity",
    "Maximum identity",
    "View details",
]

#: Fixed part of the "percent of matches" column; the operator and threshold follow it.
PERCENT_COLUMN_PREFIX = "Percent of protein sequence matches at identity"

_PERCENT_COLUMN_RE = re.compile(
    r"percent of protein sequence matches at identity"
    r"\s*(?P<op>>=|<=|>|<)?\s*(?P<threshold>[0-9]+(?:\.[0-9]+)?)?\s*%?\s*$",
    re.IGNORECASE,
)

#: Sidecar file describing how the CSVs in a directory were produced.
METADATA_FILENAME = "conservancy_meta.json"

_GREATER = {">=", ">"}
_LESSER = {"<=", "<"}


# --------------------------------------------------------------------------- threshold plumbing
def format_threshold(threshold) -> str:
    """``70.0 -> '70'``, ``72.5 -> '72.5'`` — how a threshold is spelled in a column header."""
    value = float(threshold)
    return str(int(value)) if value.is_integer() else f"{value:g}"


def percent_column(operator: str = "<=", threshold=100) -> str:
    """The "percent of matches" column header for a given operator and threshold."""
    return f"{PERCENT_COLUMN_PREFIX} {operator} {format_threshold(threshold)}%"


def operator_family(operator) -> str:
    """
    ``'ge'`` for ``>=``/``>``, ``'lt'`` for ``<``/``<=``, ``''`` when unknown.

    IEDB offers ``>=``/``<=`` while POA2 writes ``>=``/``<``; comparing families instead of the
    literal symbols keeps a ``<=`` CSV compatible with a ``-l`` (unique) run.
    """
    op = str(operator).strip()
    if op in _GREATER:
        return "ge"
    if op in _LESSER:
        return "lt"
    return ""


def parse_percent_column(name) -> Optional[Tuple[str, Optional[float]]]:
    """
    ``(operator, threshold)`` declared by a "percent of matches" header; ``None`` if not one.

    The threshold is ``None`` when the header carries no number (a hand-made CSV).
    """
    match = _PERCENT_COLUMN_RE.match(str(name).strip())
    if match is None:
        return None
    threshold = match.group("threshold")
    return (match.group("op") or "", float(threshold) if threshold is not None else None)


def find_percent_column(columns) -> str:
    """
    The "percent of matches" column, whatever threshold it was generated at.

    Raises:
        KeyError: no such column (the CSVs are not conservancy results).
        ValueError: several *different* ones — CSVs generated at different thresholds were mixed
            in the same directory, which would filter each species by a different criterion.
    """
    found = [c for c in columns if parse_percent_column(c) is not None]
    unique = list(dict.fromkeys(found))
    if not unique:
        raise KeyError(
            "ERROR: no '" + PERCENT_COLUMN_PREFIX + " ...' column in the conservancy CSVs. "
            f"Columns found: {list(columns)}"
        )
    if len(unique) > 1:
        raise ValueError(
            "ERROR: the conservancy CSVs in this directory were generated at different sequence "
            f"identity thresholds: {unique}. Re-generate them all at the same threshold, or keep "
            "each threshold in its own directory."
        )
    return unique[0]


# --------------------------------------------------------------------------- sidecar metadata
@dataclass
class ConservancyMeta:
    """How the conservancy CSVs of a directory were produced."""

    threshold: float
    operator: str
    source: str = "local"          # 'local' (our reimplementation) | 'iedb' | 'manual'
    species: Dict[str, dict] = field(default_factory=dict)
    generated_at: str = ""
    version: int = 1
    files: Optional[Dict[str, str]] = None  # relative CSV path -> SHA-256; None for legacy data

    @property
    def family(self) -> str:
        return operator_family(self.operator)


def _csv_header(path) -> Optional[List[str]]:
    """Column names of a CSV, or ``None`` when it cannot be read."""
    try:
        return list(pd.read_csv(path, nrows=0).columns)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        logger.warning("%s POA2/conservancy | could not read the header of %s: %s",
                       TAG_WARN, path, exc)
        return None


def _batch_inventory(directory: str, species: Dict[str, dict]) -> Dict[str, str]:
    """
    SHA-256 of every CSV the batch being recorded produced, keyed by path relative to ``directory``.

    The caller knows what it produced: a local run names each species' CSV, and that list is the
    inventory. Only a declaration *about* a directory (a manual upload, an IEDB download — no
    species) has to look at the directory, and even then a file that carries no "percent of
    matches" column is not a conservancy result: inventorying the researcher's own ``notas.csv``
    would turn editing it into a hard stop in :func:`inspect_directory`, and therefore in steps 4,
    5 and 7, with no way back but regenerating the whole batch. Left out, it reaches POA2 through
    :func:`stale_csvs` instead, which names it and says what to do about it.
    """
    declared = list(dict.fromkeys(str(info.get("csv")) for info in species.values()
                                  if info.get("csv")))
    if declared:
        return {name: file_hash(safe_child(directory, name)) for name in declared}
    root = Path(directory)
    return {Path(p).relative_to(root).as_posix(): file_hash(p)
            for p in map_EpConservFiles(directory)
            if any(parse_percent_column(c) is not None for c in (_csv_header(p) or []))}


def write_metadata(directory: str, threshold: float, operator: str, source: str = "local",
                   species: Optional[Dict[str, dict]] = None) -> str:
    """Record the threshold/operator that produced a directory of CSVs; returns the file path."""
    species = dict(species or {})
    threshold = float(threshold)
    if not math.isfinite(threshold) or not 0 <= threshold <= 100:
        raise ValueError("Invalid conservancy threshold: expected a finite percentage in [0, 100].")
    if not operator_family(operator):
        raise ValueError(f"Invalid conservancy operator: {operator!r}")
    meta = ConservancyMeta(
        threshold=threshold,
        operator=str(operator),
        source=source,
        species=species,
        generated_at=datetime.now().isoformat(timespec="seconds"),
        version=2,
        files=_batch_inventory(directory, species),
    )
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, METADATA_FILENAME)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(asdict(meta), fh, ensure_ascii=False, indent=2)
    return path


def _read_inventory(value) -> Optional[Dict[str, str]]:
    """
    The ``files`` mapping of a sidecar read from disk, or ``None`` when it records none.

    Raises ``TypeError`` on anything else, so that a hand-edited (or foreign-version)
    ``"files": ["A.csv"]`` is rejected where :func:`read_metadata` already handles an unusable
    sidecar, instead of escaping as an ``AttributeError`` from ``meta.files.items()`` — which the
    callers of :func:`inspect_directory` do not catch, and which reaches the user as a traceback.
    """
    if value is None:
        return None
    if not isinstance(value, dict) or not all(
            isinstance(name, str) and isinstance(digest, str) for name, digest in value.items()):
        raise TypeError("'files' must map each CSV name to its SHA-256, "
                        f"not {type(value).__name__}: {value!r}")
    return dict(value)


def read_metadata(directory: str) -> Optional[ConservancyMeta]:
    """The sidecar metadata of a CSV directory, or ``None`` when absent/unreadable."""
    path = os.path.join(directory, METADATA_FILENAME)
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return ConservancyMeta(
            threshold=float(data["threshold"]),
            operator=str(data.get("operator", "")),
            source=str(data.get("source", "local")),
            species=dict(data.get("species") or {}),
            generated_at=str(data.get("generated_at", "")),
            version=int(data.get("version", 1)),
            files=_read_inventory(data.get("files")),
        )
    except (OSError, ValueError, KeyError, TypeError) as exc:
        logger.warning("%s POA2/conservancy | ignoring an unreadable %s in %s: %s",
                       TAG_WARN, METADATA_FILENAME, directory, exc)
        return None


@dataclass
class ConservancySource:
    """What is known about the CSVs in a directory, for display and for validation."""

    csv_count: int
    threshold: Optional[float] = None
    operator: str = ""
    origin: str = "unknown"        # 'meta' | 'header' | 'unknown'
    meta: Optional[ConservancyMeta] = None

    @property
    def family(self) -> str:
        return operator_family(self.operator)


def inspect_directory(directory: str) -> ConservancySource:
    """
    Read back the threshold/operator behind a directory of conservancy CSVs.

    Metadata must agree with every declared CSV header and with recorded file hashes.
    """
    files = map_EpConservFiles(directory)
    meta = read_metadata(directory)
    declared = set()
    for path in files:
        columns = _csv_header(path)
        if columns is None:
            continue
        for column in columns:
            parsed = parse_percent_column(column)
            if parsed is not None and parsed[1] is not None:
                declared.add(parsed)
                # By family, not by symbol: an IEDB export spells the same criterion '<=' where
                # POA2 writes '<', and comparing the literals would reject a directory the rest
                # of this module — check_threshold below included — declares compatible.
                if meta is not None and (
                    abs(parsed[1] - meta.threshold) > 1e-6
                    or (parsed[0] and meta.family
                        and operator_family(parsed[0]) != meta.family)
                ):
                    raise ValueError(f"CSV header threshold/operator in {path} does not match "
                                     f"{METADATA_FILENAME}: {parsed} vs {meta.operator}{meta.threshold}.")
    if meta is not None:
        if not math.isfinite(meta.threshold) or not 0 <= meta.threshold <= 100:
            raise ValueError("Invalid threshold in conservancy metadata.")
        if meta.files is not None:
            current = {Path(p).relative_to(Path(directory)).as_posix(): p for p in files}
            for name, digest in meta.files.items():
                if name not in current or file_hash(current[name]) != digest:
                    raise ValueError(f"Conservancy CSV {name} is missing or changed since metadata was recorded.")
        return ConservancySource(len(files), meta.threshold, meta.operator, "meta", meta)
    if len(declared) == 1:
        operator, threshold = declared.pop()
        return ConservancySource(len(files), threshold, operator, "header")
    return ConservancySource(len(files))


def check_threshold(directory: str, ID_threshold, symbol, *, strict: bool = False) -> Optional[str]:
    """
    Compare the requested ``-t``/``-g|-l`` with what actually produced the CSVs.

    A mismatch means POA2 is about to filter data computed under different rules and label the
    result with the rules the user *believes* were used — the failure this check exists to stop.

    Returns:
        the mismatch description, or ``None`` when they agree (or nothing can be compared).

    Raises:
        ValueError: on a mismatch against the sidecar metadata (always), or against the CSV
            headers when ``strict``.
    """
    source = inspect_directory(directory)
    if source.origin == "unknown":
        return None

    try:
        requested = float(ID_threshold)
    except (TypeError, ValueError):
        requested = None
    requested_family = operator_family(symbol)

    problems = []
    if requested is not None and source.threshold is not None and abs(source.threshold - requested) > 1e-6:
        problems.append(
            f"threshold: POA2 was given -t {format_threshold(requested)} but the CSVs were "
            f"produced at {format_threshold(source.threshold)}"
        )
    if requested_family and source.family and requested_family != source.family:
        problems.append(
            f"operator: POA2 was given '{symbol}' but the CSVs were produced with "
            f"'{source.operator}'"
        )
    if not problems:
        return None

    where = ("the recorded metadata (" + METADATA_FILENAME + ")" if source.origin == "meta"
             else "the CSV headers")
    message = (
        "ERROR: the sequence identity threshold requested does not match "
        f"{where} in '{directory}' — " + "; ".join(problems) + ". "
        "Re-run the Epitope Conservancy Analysis at the intended threshold, or run POA2 with the "
        "values the CSVs were actually produced with."
    )
    if source.origin == "meta" or strict:
        raise ValueError(message)
    logger.warning("%s", message)
    return message


# --------------------------------------------------------------------------- CSV helpers
def map_EpConservFiles(directory):
    """Retrieve file paths for all CSV files in ``directory`` and subdirectories."""
    EpConservFiles = []
    for folders, subfolders, files in os.walk(directory):
        subfolders[:] = [d for d in subfolders if not d.startswith(".poa_")]
        for file in files:
            path_file = os.path.join(folders, file)
            path_file = str(path_file)
            if file.lower().endswith(".csv"):
                EpConservFiles.append(path_file)
    return sorted(EpConservFiles)


def effective_match_threshold(seq_match: float, n_proteins: int) -> Optional[float]:
    """
    The value ``-m`` actually behaves as, given a comparison set of ``n_proteins``.

    The percent column is ``matched/n``, so it can only ever take the ``n + 1`` values
    ``0, 100/n, 200/n, ..., 100``. A ``-m`` that falls between two of them behaves as the next
    one up — with ``n = 2`` the reachable values are 0, 50 and 100, so the default ``-m 60``
    is indistinguishable from demanding **every** protein match. The researcher picked a number
    meaning "most of them" and got unanimity.

    Args:
        seq_match: the ``-m`` the user chose, in percent.
        n_proteins: size of the comparison set.

    Returns:
        The smallest reachable percentage that is ``>= seq_match``, or ``None`` when ``n`` is
        unknown/zero (nothing can be said).
    """
    if not n_proteins or n_proteins <= 0:
        return None
    step = 100.0 / n_proteins
    for k in range(n_proteins + 1):
        value = k * step
        if value >= seq_match - 1e-9:
            return round(value, 6)
    return 100.0


def check_match_resolution(directory: str, seq_match: float) -> List[str]:
    """
    Report species whose comparison set is too small for the chosen ``-m`` to mean what it says.

    Complements the ``len(pairs) < 2`` warning raised when the CSVs are produced: that one only
    fires for a single protein, but the gap between the chosen ``-m`` and its effective value
    stays wide for several proteins (``n = 3`` still jumps 0 → 33 → 67 → 100).

    Args:
        directory: the folder of conservancy CSVs.
        seq_match: POA2's ``-m``.

    Returns:
        One human-readable line per affected species, empty when every set is fine or there is
        no sidecar to tell how large they were.
    """
    meta = read_metadata(directory)
    if meta is None or not meta.species:
        return []

    avisos = []
    for specie, info in sorted(meta.species.items()):
        n = int(info.get("n_proteins") or 0)
        effective = effective_match_threshold(seq_match, n)
        if effective is None or effective <= seq_match + 1e-9:
            continue
        detalhe = ("equivale a exigir que TODAS correspondam"
                   if effective >= 100.0 else f"equivale a -m {effective:g}")
        avisos.append(
            f"{specie}: o conjunto tem {n} proteína(s), então a conservância só assume valores "
            f"de {100.0 / n:g} em {100.0 / n:g} pontos — o -m {seq_match:g} escolhido "
            f"{detalhe}.")
    return avisos


def stale_csvs(directory: str) -> List[str]:
    """
    CSVs in ``directory`` that the run described by the sidecar did not produce.

    Nothing ever empties the conservancy directory, and :func:`EpitConservAnalysis` concatenates
    **every** CSV it finds under it. So a species dropped between two runs — or a leftover from
    an earlier analysis — keeps being folded into the result, and when the header happens to
    match it is completely invisible.

    Args:
        directory: the ``-d`` of POA2.

    Returns:
        File names (not paths) present on disk but absent from the sidecar's species list,
        sorted. Empty when there is no sidecar, since then there is nothing to compare against
        and :func:`check_threshold` already reports the missing provenance.
    """
    meta = read_metadata(directory)
    if meta is None:
        return []
    if meta.files is not None:
        produced = set(meta.files)
    elif meta.species:
        produced = {str(info.get("csv", "")) for info in meta.species.values()}
    else:
        return []  # legacy manual metadata had no inventory, not an empty inventory
    found = {Path(p).relative_to(Path(directory)).as_posix() for p in map_EpConservFiles(directory)}
    return sorted(found - produced)


def strip_percent(series: pd.Series) -> pd.Series:
    """``'80.00%' -> 80.0``; non-numeric and missing cells become NaN instead of raising."""
    text = series.astype(str).str.strip().str.rstrip("%").str.strip()
    return pd.to_numeric(text, errors="coerce")


# --------------------------------------------------------------------------- main filter
def EpitConservAnalysis(ID_threshold, symbol, seq_match, max_ID, min_ID, directory,
                        *, identity_filter: bool = False, strict: bool = False):
    """
    Processes and filters the results of Epitope Conservancy Analysis from multiple CSV files.

    Parameters:
        ID_threshold: sequence-identity threshold the CSVs were produced with (POA2 ``-t``).
        symbol: ``'>='`` (conserved) or ``'<'`` (unique) — the operator that goes with it.
        seq_match: minimum percent of protein matches to keep an epitope (POA2 ``-m``).
        max_ID / min_ID: bounds on each epitope's maximum/minimum identity.
        directory: folder holding the conservancy CSVs.
        identity_filter: additionally require the epitope's own identity to satisfy
            ``{symbol} ID_threshold``. Off by default because it changes which epitopes POA2
            selects: ``>=`` requires the *minimum* identity to reach the threshold (conserved
            across the whole set), ``<`` requires the *maximum* to stay below it (unique).
        strict: raise, rather than warn, when the CSV headers declare a different threshold.

    Returns:
        pd.DataFrame: filtered and processed conservancy analysis results.
    """
    # Refuse to filter data that was not produced under the requested threshold/operator.
    check_threshold(directory, ID_threshold, symbol, strict=strict)

    # Every CSV under the directory is concatenated below, so one left over from an earlier run
    # silently joins this analysis. Under -strict that is a stop; otherwise it is named in the
    # log, because deleting a file the user put there is not this function's call to make.
    leftovers = stale_csvs(directory)
    if leftovers:
        message = (f"{len(leftovers)} CSV(s) in {directory} were not produced by the run recorded "
                   f"in {METADATA_FILENAME} and will be merged into this analysis anyway: "
                   f"{', '.join(leftovers)}. Remove them, or re-run the conservancy step.")
        if strict:
            raise ValueError(message)
        logger.warning(message)

    # Not an error — a small comparison set is a fact about the data, not a mistake in the run —
    # but the gap between the -m that was typed and the one that takes effect has to be visible.
    for aviso in check_match_resolution(directory, seq_match):
        logger.warning("%s", aviso)

    # Get the list of CSV files in the directory
    filesList = map_EpConservFiles(directory)
    if not filesList:
        return pd.DataFrame(columns=IEDB_COLUMNS)

    # Consolidate data from all CSV files. (Was an outer pd.merge against an empty seed frame with
    # hardcoded column names, which duplicated the percent column — leaving it all-NaN — as soon as
    # the CSVs carried a threshold other than the default '<= 100%'.) The merge also collapsed
    # byte-identical rows, which only happen when the same CSV is present twice; drop_duplicates
    # keeps that, so a duplicated file still does not double an epitope.
    allSp_ConservDF = pd.concat([pd.read_csv(file) for file in filesList], ignore_index=True)
    allSp_ConservDF = allSp_ConservDF.drop_duplicates(ignore_index=True)

    percent_col = find_percent_column(allSp_ConservDF.columns)

    # Clean and preprocess the data
    allSp_ConservDF[percent_col] = (
        allSp_ConservDF[percent_col].astype(str).str.replace("%", "", regex=False)
    )
    allSp_ConservDF["Minimum identity"] = strip_percent(allSp_ConservDF["Minimum identity"])
    allSp_ConservDF["Maximum identity"] = strip_percent(allSp_ConservDF["Maximum identity"])

    # Apply filters based on minimum and maximum identity thresholds
    allSp_ConservDF_idMin = allSp_ConservDF[allSp_ConservDF["Minimum identity"] >= min_ID]
    if allSp_ConservDF_idMin.empty:  # Check if the DataFrame is empty after filtering
        return allSp_ConservDF_idMin

    allSp_ConservDF_idMax = allSp_ConservDF_idMin[allSp_ConservDF_idMin["Maximum identity"] <= max_ID].copy()
    if allSp_ConservDF_idMax.empty:  # Check if the DataFrame is empty after filtering
        return allSp_ConservDF_idMax

    # Split and process the "Percent of protein sequence matches" column
    seqMatchSplit = allSp_ConservDF_idMax[percent_col].str.split(" ", n=1, expand=True)
    allSp_ConservDF_idMax["Percent"] = pd.to_numeric(seqMatchSplit[0], errors="coerce")
    allSp_ConservDF_idMax["Reason"] = seqMatchSplit[1] if seqMatchSplit.shape[1] > 1 else ""

    # Filter based on the minimum sequence match percentage
    allSp_ConservDF_final = allSp_ConservDF_idMax[allSp_ConservDF_idMax["Percent"] >= seq_match]

    # Optional: make the threshold a selection criterion, not just the basis of the percent column
    if identity_filter:
        allSp_ConservDF_final = _apply_identity_filter(allSp_ConservDF_final, ID_threshold, symbol)

    # Rename columns for clarity
    new_colNames = {
        percent_col: f"Protein(s) sequence match(es) at Sequence identity threshold {symbol}{ID_threshold}",
        "Minimum identity": "Minimum identity(%)",
        "Maximum identity": "Maximum identity(%)",
    }
    allSp_ConservDF_final = allSp_ConservDF_final.rename(columns=new_colNames, inplace=False)

    # Select and return the final columns of interest
    final_Conserv_df = allSp_ConservDF_final[[
        "Epitope name", "Epitope sequence", "Epitope length", f"Protein(s) sequence match(es) at Sequence identity threshold {symbol}{ID_threshold}", "Minimum identity(%)", "Maximum identity(%)"
    ]]

    return final_Conserv_df


def stage_counts(directory, ID_threshold=70, symbol=">=", seq_match=0, max_ID=100, min_ID=0,
                 identity_filter=False) -> List[Dict]:
    """
    How many epitopes survive each filter of :func:`EpitConservAnalysis`, in order.

    POA2 reports only what came through every filter, so a run that selects nothing looks the
    same whether the identity bounds were too tight or no epitope reached ``-m``. This reapplies
    the same steps cumulatively and counts what is left after each one.

    This mirrors the filter chain in :func:`EpitConservAnalysis` and the two must stay in step —
    ``tests/test_step_outputs.py`` asserts the final count equals ``len(EpitConservAnalysis(...))``
    for the same arguments, so a change to one that is not made to the other fails there.

    Returns one dict per stage — ``{'stage', 'value', 'remaining'}``, where ``stage`` is a
    language-neutral key for the caller to label — or ``[]`` when the directory holds no CSV, or
    CSVs this cannot read. Unlike :func:`EpitConservAnalysis` it neither raises on a threshold
    mismatch nor on an unreadable header: it describes the data rather than deciding anything
    with it, and a description that cannot be produced must not stop the run it describes.
    """
    files = map_EpConservFiles(directory)
    if not files:
        return []

    try:
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        df = df.drop_duplicates(ignore_index=True)
        percent_col = find_percent_column(df.columns)
    except (OSError, ValueError, KeyError, pd.errors.ParserError):
        return []
    if not {"Minimum identity", "Maximum identity"} <= set(df.columns):
        return []
    df["Minimum identity"] = strip_percent(df["Minimum identity"])
    df["Maximum identity"] = strip_percent(df["Maximum identity"])

    stages = [{"stage": "input", "value": None, "remaining": len(df)}]

    df = df[df["Minimum identity"] >= min_ID]
    stages.append({"stage": "min_identity", "value": min_ID, "remaining": len(df)})

    df = df[df["Maximum identity"] <= max_ID].copy()
    stages.append({"stage": "max_identity", "value": max_ID, "remaining": len(df)})

    # The percent column reads '100.00% (1/1)'; only the leading number is the match percentage.
    if not df.empty:
        percent = df[percent_col].astype(str).str.replace("%", "", regex=False)
        df["Percent"] = pd.to_numeric(percent.str.split(" ", n=1, expand=True)[0], errors="coerce")
        df = df[df["Percent"] >= seq_match]
    stages.append({"stage": "seq_match", "value": seq_match, "remaining": len(df)})

    if identity_filter:
        if not df.empty:
            df = _apply_identity_filter(df, ID_threshold, symbol)
        stages.append({"stage": "identity_filter", "value": (symbol, ID_threshold),
                       "remaining": len(df)})

    return stages


def _apply_identity_filter(dataframe: pd.DataFrame, ID_threshold, symbol) -> pd.DataFrame:
    """
    Keep only epitopes whose own identity satisfies ``{symbol} ID_threshold``.

    Conserved (``>=``) looks at the **minimum** identity: the epitope must reach the threshold
    against every protein of the set. Unique (``<``) looks at the **maximum**: it must stay below
    the threshold against all of them.
    """
    family = operator_family(symbol)
    try:
        threshold = float(ID_threshold)
    except (TypeError, ValueError):
        logger.warning("%s POA2/conservancy | identity filter skipped: '%s' is not a usable "
                       "threshold.", TAG_WARN, ID_threshold)
        return dataframe
    if family == "ge":
        return dataframe[dataframe["Minimum identity"] >= threshold]
    if family == "lt":
        return dataframe[dataframe["Maximum identity"] < threshold]
    logger.warning("%s POA2/conservancy | identity filter skipped: unknown operator '%s'.",
                   TAG_WARN, symbol)
    return dataframe
