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
import os
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Dict, Optional, Tuple

import pandas as pd

from ...logging_conf import get_logger

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

    @property
    def family(self) -> str:
        return operator_family(self.operator)


def write_metadata(directory: str, threshold, operator: str, source: str = "local",
                   species: Optional[Dict[str, dict]] = None) -> str:
    """Record the threshold/operator that produced a directory of CSVs; returns the file path."""
    meta = ConservancyMeta(
        threshold=float(threshold),
        operator=str(operator),
        source=source,
        species=dict(species or {}),
        generated_at=datetime.now().isoformat(timespec="seconds"),
    )
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, METADATA_FILENAME)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(asdict(meta), fh, ensure_ascii=False, indent=2)
    return path


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
        )
    except (OSError, ValueError, KeyError, TypeError) as exc:
        logger.warning("Ignoring unreadable %s in %s: %s", METADATA_FILENAME, directory, exc)
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

    The sidecar metadata wins when present (we wrote it, so it is unambiguous); otherwise the
    header of the CSVs themselves is parsed, which is what an IEDB download carries.
    """
    files = map_EpConservFiles(directory)
    meta = read_metadata(directory)
    if meta is not None:
        return ConservancySource(len(files), meta.threshold, meta.operator, "meta", meta)

    declared = set()
    for path in files:
        try:
            columns = pd.read_csv(path, nrows=0).columns
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            logger.warning("Could not read the header of %s: %s", path, exc)
            continue
        for column in columns:
            parsed = parse_percent_column(column)
            if parsed is not None and parsed[1] is not None:
                declared.add(parsed)
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
        for file in files:
            path_file = os.path.join(folders, file)
            path_file = str(path_file)
            if ".csv" in file:
                EpConservFiles.append(path_file)
    return EpConservFiles


def _strip_percent(series: pd.Series) -> pd.Series:
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
    allSp_ConservDF["Minimum identity"] = _strip_percent(allSp_ConservDF["Minimum identity"])
    allSp_ConservDF["Maximum identity"] = _strip_percent(allSp_ConservDF["Maximum identity"])

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
        logger.warning("Identity filter skipped: '%s' is not a usable threshold.", ID_threshold)
        return dataframe
    if family == "ge":
        return dataframe[dataframe["Minimum identity"] >= threshold]
    if family == "lt":
        return dataframe[dataframe["Maximum identity"] < threshold]
    logger.warning("Identity filter skipped: unknown operator '%s'.", symbol)
    return dataframe
