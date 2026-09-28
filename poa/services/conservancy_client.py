"""Local reimplementation of the IEDB Epitope Conservancy Analysis.

IEDB provides NO API and NO downloadable standalone for this tool (verified 2026), so — with the
user's approval — this reproduces the algorithm locally (Bui et al., BMC Bioinformatics 2007):

    For a linear epitope of length L, slide an ungapped window of length L across each protein
    sequence and compute the percent identity (matches / L * 100) at every offset; the epitope's
    identity to that protein is the MAXIMUM over all offsets. Given a sequence-identity threshold
    T and an operator, the tool reports, per epitope:
      * the percent (and count) of proteins that contain a match with identity {>= | <} T,
      * the minimum and maximum per-protein identity across the whole protein set.

The output CSV uses the column headers that :mod:`poa.core.parsers.conservancy` consumes, with the
threshold and operator spelled out in the "percent of matches" header exactly as IEDB writes them,
so POA2 works identically whether the CSVs come from IEDB or from here — and can check that it is
filtering data produced under the rules it was told about. A ``conservancy_meta.json`` sidecar
records the same thing unambiguously.

IMPORTANT: validate this against a real IEDB CSV export for your dataset before relying on it in
production — see ``tests/test_conservancy_local.py`` for the equivalence checks.
"""
from __future__ import annotations

import csv
import io
import os
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from Bio import SeqIO

from ..core.parsers.conservancy import operator_family, percent_column, write_metadata
from ..logging_conf import get_logger

logger = get_logger("conservancy")

#: Operators the analysis accepts: conserved (``>=``) and unique (``<``).
OPERATORS = (">=", "<")


def csv_columns(operator: str = "<=", threshold=100) -> List[str]:
    """Headers of a conservancy CSV produced at ``operator``/``threshold``."""
    return [
        "Epitope #",
        "Epitope name",
        "Epitope sequence",
        "Epitope length",
        percent_column(operator, threshold),
        "Minimum identity",
        "Maximum identity",
        "View details",
    ]


#: Legacy alias — the IEDB default header (``<= 100%``), kept for callers that predate the
#: threshold-aware header.
CSV_COLUMNS = csv_columns()


def max_identity(epitope: str, protein: str) -> float:
    """Maximum ungapped percent identity of ``epitope`` against any window of ``protein``."""
    epitope = epitope.upper()
    protein = protein.upper()
    L = len(epitope)
    if L == 0 or len(protein) < L:
        return 0.0
    best = 0.0
    for offset in range(len(protein) - L + 1):
        window = protein[offset:offset + L]
        matches = sum(1 for a, b in zip(epitope, window) if a == b)
        identity = matches / L * 100.0
        if identity > best:
            best = identity
            if best == 100.0:
                break
    return best


def matches_threshold(identity: float, threshold: float, operator: str = ">=") -> bool:
    """Whether a per-protein identity counts as a match under ``operator``/``threshold``."""
    family = operator_family(operator)
    if family == "lt":
        return identity < threshold
    if family == "ge":
        return identity >= threshold
    raise ValueError(f"Unsupported conservancy operator: {operator!r} (use one of {OPERATORS}).")


def conservancy_rows(
    epitopes: Sequence[Tuple[str, str]],
    proteins: Sequence[Tuple[str, str]],
    threshold: float,
    operator: str = ">=",
) -> List[dict]:
    """
    Compute conservancy for each epitope against the protein set at ``threshold`` percent identity.

    Parameters:
        epitopes: list of (name, sequence).
        proteins: list of (name, sequence).
        threshold: sequence-identity threshold (percent) used for the "matches" count.
        operator: ``'>='`` counts proteins the epitope is conserved in, ``'<'`` counts those it
            stays below — the "unique epitopes" objective, which previously had no effect here
            because the comparison was hardcoded to ``>=``.

    Returns:
        list of dict rows keyed by ``csv_columns(operator, threshold)``.
    """
    percent_col = percent_column(operator, threshold)
    n_proteins = len(proteins)
    rows = []
    for i, (name, seq) in enumerate(epitopes, start=1):
        per_protein = [max_identity(seq, pseq) for _, pseq in proteins]
        matched = sum(1 for v in per_protein if matches_threshold(v, threshold, operator))
        percent = (matched / n_proteins * 100.0) if n_proteins else 0.0
        min_id = min(per_protein) if per_protein else 0.0
        max_id = max(per_protein) if per_protein else 0.0
        rows.append({
            "Epitope #": i,
            "Epitope name": name,
            "Epitope sequence": seq,
            "Epitope length": len(seq),
            percent_col: f"{percent:.2f}% ({matched}/{n_proteins})",
            "Minimum identity": f"{min_id:.2f}%",
            "Maximum identity": f"{max_id:.2f}%",
            "View details": "-",
        })
    return rows


def rows_to_csv(rows: List[dict], columns: Optional[Sequence[str]] = None) -> str:
    """Serialize conservancy rows to an IEDB-compatible CSV string."""
    fieldnames = list(columns) if columns else (list(rows[0].keys()) if rows else CSV_COLUMNS)
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return buf.getvalue()


def _read_fasta_pairs(fasta_path: str) -> List[Tuple[str, str]]:
    return [(rec.id, str(rec.seq)) for rec in SeqIO.parse(fasta_path, "fasta")]


def run_conservancy(
    epitopes_fasta: str,
    proteins: str | Sequence[Tuple[str, str]],
    threshold: float,
    operator: str = ">=",
) -> str:
    """
    Run the local conservancy analysis for a single epitope FASTA against a protein set.

    Parameters:
        epitopes_fasta: path to a FASTA of epitopes (e.g. a POA1 ``<sp>_epitopes.fasta``).
        proteins: path to a protein FASTA, or a list of (name, sequence).
        threshold: sequence-identity threshold (percent).
        operator: ``'>='`` (conserved) or ``'<'`` (unique).

    Returns:
        CSV text (IEDB-compatible).
    """
    epitopes = _read_fasta_pairs(epitopes_fasta)
    protein_pairs = _read_fasta_pairs(proteins) if isinstance(proteins, str) else list(proteins)
    rows = conservancy_rows(epitopes, protein_pairs, threshold, operator)
    return rows_to_csv(rows, csv_columns(operator, threshold))


def specie_of_header(name: str) -> str:
    """Species token of a ``Protein_Specie_ID`` header (upper-cased); ``''`` if absent."""
    parts = str(name).upper().split("_")
    return parts[1] if len(parts) > 1 else ""


def group_proteins_by_specie(
    proteins: Sequence[Tuple[str, str]]
) -> Dict[str, List[Tuple[str, str]]]:
    """Bucket ``(name, sequence)`` pairs by the species token of their header."""
    groups: Dict[str, List[Tuple[str, str]]] = {}
    for name, seq in proteins:
        groups.setdefault(specie_of_header(name), []).append((name, seq))
    groups.pop("", None)
    return groups


def run_conservancy_for_dir(
    poa1_conservancy_dir: str,
    proteins: str | Sequence[Tuple[str, str]],
    threshold: float,
    out_dir: str,
    specie_proteins: Optional[Mapping[str, "str | Sequence[Tuple[str, str]]"]] = None,
    per_specie: bool = True,
    operator: str = ">=",
) -> List[str]:
    """
    Process every ``*_epitopes.fasta`` produced by POA1 and write one conservancy CSV per file.

    This replaces the manual "submit each FASTA to the IEDB web tool and download the CSV" step.

    POA1 writes one ``<SPECIE>_epitopes.fasta`` per species, and each species must be compared
    against **its own** protein set — comparing DENV1 epitopes to a CHIKV protein set yields
    meaningless identities. The comparison set for ``<SPECIE>`` is chosen in this order:

    1. ``specie_proteins[SPECIE]`` — an explicit per-species set (e.g. its world/diversity FASTA);
    2. the records of ``proteins`` whose header species token is ``SPECIE``;
    3. the whole ``proteins`` set (with a warning when it covers other species too).

    A ``conservancy_meta.json`` is written alongside the CSVs recording the threshold, the operator
    and, per species, the size of the comparison set — which POA2 reads back to confirm it is
    filtering what it thinks it is filtering.

    Returns:
        list of written CSV file paths (ready to be consumed by POA2 via ``-d``).
    """
    protein_pairs = _read_fasta_pairs(proteins) if isinstance(proteins, str) else list(proteins)
    by_specie = group_proteins_by_specie(protein_pairs) if per_specie else {}
    explicit = {
        str(k).upper(): (_read_fasta_pairs(v) if isinstance(v, str) else list(v))
        for k, v in (specie_proteins or {}).items()
    }
    os.makedirs(out_dir, exist_ok=True)

    written = []
    species_meta = {}
    for fname in sorted(os.listdir(poa1_conservancy_dir)):
        if not fname.endswith("_epitopes.fasta"):
            continue
        specie = fname[: -len("_epitopes.fasta")].upper()

        if explicit.get(specie):
            pairs, origin = explicit[specie], "conjunto próprio da espécie"
        elif by_specie.get(specie):
            pairs, origin = by_specie[specie], f"{len(by_specie[specie])} proteína(s) de {specie}"
        else:
            pairs, origin = protein_pairs, "conjunto completo (sem correspondência por espécie)"
            if len(by_specie) > 1:
                logger.warning(
                    "%s: no protein of species '%s' in the comparison set (species found: %s). "
                    "Falling back to the full set — the identities will mix species.",
                    fname, specie, ", ".join(sorted(by_specie)),
                )

        # With a single protein to compare against, every epitope matches it (the epitope was
        # predicted ON that protein), so the percent column is 100% whatever the threshold is.
        # Silently reporting "100% conserved" is the failure this warning exists to make visible.
        if len(pairs) < 2:
            logger.warning(
                "%s: the comparison set for '%s' has %d protein(s), so the %s%g%% threshold cannot "
                "discriminate — every epitope will come out at 100%%. Supply a diversity set "
                "(e.g. the species' world set) for a meaningful conservancy.",
                fname, specie, len(pairs), operator, threshold,
            )

        src = os.path.join(poa1_conservancy_dir, fname)
        epitopes = _read_fasta_pairs(src)
        rows = conservancy_rows(epitopes, pairs, threshold, operator)
        csv_text = rows_to_csv(rows, csv_columns(operator, threshold))
        out_name = fname.replace("_epitopes.fasta", "_conservancy.csv")
        out_path = os.path.join(out_dir, out_name)
        with open(out_path, "w", encoding="utf-8", newline="") as fh:
            fh.write(csv_text)
        logger.info("Conservancy CSV written: %s (%s vs %s)", out_path, specie, origin)
        written.append(out_path)
        species_meta[specie] = {
            "csv": out_name,
            "n_proteins": len(pairs),
            "n_epitopes": len(epitopes),
            "origin": origin,
            "discriminating": len(pairs) > 1,
        }

    if written:
        write_metadata(out_dir, threshold, operator, source="local", species=species_meta)
    return written
