"""IEDB MHC-II Binding Predictions client (legacy REST API, synchronous).

Endpoint (verified live, 2026): https://tools-cluster-interface.iedb.org/tools_api/mhcii/
POST form fields: method, sequence_text, allele (comma-separated), length. Returns a TSV table.

The response is normalized and written as ``<Protein>_<Specie>.html`` so the existing
:mod:`poa.core.parsers.mhcii` parser consumes it unchanged (its "manual upload" path uses the
very same parser). On any network/HTTP problem a :class:`ServiceUnavailable` is raised so the
caller can fall back to a manually downloaded result file.
"""
from __future__ import annotations

import csv
import io
import os
from typing import List, Optional, Sequence

from ..logging_conf import get_logger, log_step
from .base import Cache, ServiceError, ServiceResult, ServiceUnavailable, cached_call

logger = get_logger("mhcii")

IEDB_MHCII_API = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# Columns required by poa.core.parsers.mhcii.MHCIIAntigenEpitopes
REQUIRED_COLUMNS = [
    "allele", "seq_num", "start", "end", "method", "peptide",
    "smm_align_ic50", "nn_align_ic50", "nn_align_rank", "nn_align_adjusted_rank",
]

# Best-effort mapping of alternative IEDB column names -> the names the parser expects.
_ALIASES = {
    "nn_align_ic50": ["nn_align_ic50", "ic50", "nn_ic50"],
    "nn_align_rank": ["nn_align_rank", "rank", "percentile_rank", "nn_rank"],
    "nn_align_adjusted_rank": ["nn_align_adjusted_rank", "adjusted_rank", "nn_adjusted_rank"],
    "smm_align_ic50": ["smm_align_ic50", "smm_ic50"],
    "method": ["method"],
    "seq_num": ["seq_num", "seq_no", "sequence_number"],
    "start": ["start"],
    "end": ["end"],
    "peptide": ["peptide"],
    "allele": ["allele"],
}


def submit(sequence_text: str, alleles: Sequence[str] | str, length=15,
           method: str = "nn_align", timeout: int = 180) -> str:
    """POST to the IEDB MHC-II API and return the raw TSV text. Raises ServiceUnavailable on failure."""
    try:
        import requests
    except ImportError as exc:  # pragma: no cover - requests is a declared dependency
        raise ServiceUnavailable("The 'requests' package is required for the IEDB MHC-II API.") from exc

    allele_str = ",".join(alleles) if not isinstance(alleles, str) else alleles
    data = {"method": method, "sequence_text": sequence_text, "allele": allele_str, "length": str(length)}
    try:
        resp = requests.post(IEDB_MHCII_API, data=data, timeout=timeout)
    except requests.RequestException as exc:
        raise ServiceUnavailable(f"Could not reach the IEDB MHC-II API: {exc}") from exc

    if resp.status_code != 200:
        raise ServiceUnavailable(f"IEDB MHC-II API returned HTTP {resp.status_code}: {resp.text[:200]}")

    text = resp.text
    lowered = text.lstrip().lower()
    if not text.strip() or lowered.startswith("<html") or lowered.startswith("<!doctype"):
        raise ServiceError("Unexpected (non-tabular) response from the IEDB MHC-II API.")
    if "allele" not in text.splitlines()[0].lower():
        raise ServiceError("IEDB MHC-II response is missing the expected 'allele' header.")
    return text


def _find_header_index(lines: List[str]) -> int:
    for i, line in enumerate(lines):
        if line.lower().startswith("allele"):
            return i
    raise ServiceError("No 'allele' header line found in the MHC-II result.")


def normalize_tsv(tsv_text: str, requested_method: str = "nn_align") -> str:
    """
    Rewrite an IEDB MHC-II TSV so it contains exactly the columns the POA parser needs
    (aliasing known column variants; missing columns filled with '-').
    """
    lines = [ln for ln in tsv_text.splitlines() if ln.strip()]
    hidx = _find_header_index(lines)
    header = lines[hidx].split("\t")
    reader = csv.DictReader(io.StringIO("\n".join(lines[hidx:])), delimiter="\t")

    present = {h.lower(): h for h in header}

    def resolve(col: str) -> Optional[str]:
        for alias in _ALIASES.get(col, [col]):
            if alias.lower() in present:
                return present[alias.lower()]
        return None

    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(REQUIRED_COLUMNS)
    for row in reader:
        new_row = []
        for col in REQUIRED_COLUMNS:
            src = resolve(col)
            value = row.get(src, "-") if src else "-"
            if col == "method" and (value in (None, "-", "")):
                value = requested_method
            new_row.append(value if value not in (None, "") else "-")
        writer.writerow(new_row)
    return out.getvalue()


def write_result_file(tsv_text: str, protein: str, specie: str, out_dir: str,
                      normalize: bool = True, requested_method: str = "nn_align") -> str:
    """Write the (normalized) TSV as ``<Protein>_<Specie>.html`` for the directory-based parser."""
    os.makedirs(out_dir, exist_ok=True)
    content = normalize_tsv(tsv_text, requested_method) if normalize else tsv_text
    path = os.path.join(out_dir, f"{protein}_{specie}.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return path


def predict(sequence_text: str, alleles: Sequence[str] | str, length=15, method: str = "nn_align",
            cache: Optional[Cache] = None, use_cache: bool = True, timeout: int = 180) -> ServiceResult:
    """Cached MHC-II prediction. Returns a :class:`ServiceResult` with raw TSV content."""
    allele_str = ",".join(alleles) if not isinstance(alleles, str) else alleles
    params = {"method": method, "allele": allele_str, "length": length}
    with log_step(logger, "service/iedb-mhcii") as step:
        step.note("method=%s alleles=%s length=%s", method, allele_str, length)
        result = cached_call(
            cache, "mhcii", "iedb-api", params, [sequence_text],
            lambda: submit(sequence_text, allele_str, length, method, timeout),
            source="api", ext="tsv", use_cache=use_cache,
        )
        step.result(source=result.source, rows=max(len((result.content or "").splitlines()) - 1, 0))
    return result
