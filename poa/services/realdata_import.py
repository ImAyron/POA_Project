"""Adapter for the original manual-workflow research files.

The real research inputs (BepiPred-2.0 JSON exports, world/diversity FASTA sets) use file and
header conventions that predate the ``poa`` package and do **not** match what the pipeline
parsers expect. This module adapts the *naming* only — it never changes the predicted data:

* BepiPred-2.0 JSON exports carry a single generic antigen key (``"Sequence"``). The
  :func:`poa.core.parsers.bepipred.bp2_AntigenEpitopes` parser extracts Protein/Specie/ID from
  that key via ``(\\w+?)_(\\w+?)_(\\w+)``, so a generic key yields ``NaN`` species. We rewrite the
  key to a conforming ``Protein_Specie_ID`` string.
* POA1 needs a protein FASTA (``-f``). The exact sequence BepiPred ran on is embedded in the
  JSON's own ``AA`` array, so we rebuild the reference FASTA from it (same header convention).

The world FASTA sets need no adaptation: :mod:`poa.services.conservancy_client` reads them
directly and tolerates ambiguous ``X`` residues (they simply fail to match, lowering identity).

Note: :func:`poa.core.report.checkIntegrity` rejects any ``-f`` containing ``X``. Some references
themselves contain ``X`` (e.g. DENV3 in the sample data); such a virus is reported with
``has_x=True`` and will be refused by POA1 until the ambiguous residue is resolved upstream.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from ..logging_conf import get_logger

logger = get_logger("realdata_import")

# bepipred_<virus>.json  ->  <virus>
_FNAME_RE = re.compile(r"bepipred[_-]?(?P<virus>[a-zA-Z0-9]+)", re.IGNORECASE)


@dataclass
class Prepared:
    """A single virus adapted to pipeline conventions."""

    specie: str
    protein: str
    id_seq: str
    sequence: str
    reference_fasta: str  # path to the rebuilt -f FASTA
    b2_json: str          # path to the rewritten -b2 JSON
    has_x: bool           # reference contains an ambiguous 'X' (POA1 will refuse it)

    @property
    def header(self) -> str:
        return f"{self.protein}_{self.specie}_{self.id_seq}"


def specie_from_filename(path: str) -> str:
    """Guess the species token from a ``bepipred_<virus>.json`` filename (upper-cased)."""
    stem = os.path.splitext(os.path.basename(path))[0]
    m = _FNAME_RE.search(stem)
    return (m.group("virus") if m else stem).upper()


def _antigen_sequence(antigen: dict) -> str:
    """Join the per-residue ``AA`` array back into the reference sequence string."""
    return "".join(antigen["AA"])


def prepare_bepipred2(
    json_path: str,
    specie: str,
    out_dir: str,
    protein: str = "E",
    id_seq: str = "ref",
) -> Prepared:
    """
    Adapt one BepiPred-2.0 JSON export to pipeline conventions.

    Writes, into ``out_dir``:
      * ``<specie>_ref.fasta`` — the reference protein (``-f``), header ``Protein_Specie_ID``.
      * ``<specie>_b2.json``   — a copy whose antigen key(s) follow ``Protein_Specie_ID`` so the
                                 BepiPred-2.0 parser extracts species/protein correctly.

    The predicted per-residue scores are copied verbatim; only the antigen key changes.
    """
    with open(json_path, "r") as fh:
        data = json.load(fh)

    antigens = data.get("antigens", {})
    if not antigens:
        raise ValueError(f"{json_path}: no 'antigens' block found")

    os.makedirs(out_dir, exist_ok=True)
    items = list(antigens.items())
    multi = len(items) > 1

    new_antigens: Dict[str, dict] = {}
    records: List[Tuple[str, str]] = []
    for i, (_old_key, antigen) in enumerate(items, start=1):
        this_id = f"{id_seq}{i}" if multi else id_seq
        new_key = f"{protein}_{specie}_{this_id}"
        new_antigens[new_key] = antigen
        records.append((new_key, _antigen_sequence(antigen)))
    data["antigens"] = new_antigens

    b2_path = os.path.join(out_dir, f"{specie}_b2.json")
    with open(b2_path, "w") as fh:
        json.dump(data, fh)

    ref_path = os.path.join(out_dir, f"{specie}_ref.fasta")
    with open(ref_path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")

    full_seq = "".join(seq for _, seq in records)
    has_x = "x" in full_seq.lower()
    if has_x:
        logger.warning("%s: reference for %s contains 'X' — POA1 will refuse it (-f integrity check).", json_path, specie)

    return Prepared(
        specie=specie,
        protein=protein,
        id_seq=id_seq,
        sequence=full_seq,
        reference_fasta=ref_path,
        b2_json=b2_path,
        has_x=has_x,
    )


def prepare_bepipred2_dir(
    pred_dir: str,
    out_dir: str,
    protein: str = "E",
    protein_map: Optional[Dict[str, str]] = None,
) -> List[Prepared]:
    """
    Adapt every ``bepipred_*.json`` in ``pred_dir``.

    Parameters:
        protein: default protein label for the ``Protein_Specie_ID`` header.
        protein_map: optional per-species override, keyed by the upper-cased species token
                     (e.g. ``{"CHIKV": "E1"}``).
    """
    protein_map = {k.upper(): v for k, v in (protein_map or {}).items()}
    prepared: List[Prepared] = []
    for fname in sorted(os.listdir(pred_dir)):
        if not (fname.lower().startswith("bepipred") and fname.lower().endswith(".json")):
            continue
        json_path = os.path.join(pred_dir, fname)
        specie = specie_from_filename(json_path)
        prot = protein_map.get(specie, protein)
        prepared.append(prepare_bepipred2(json_path, specie, out_dir, protein=prot))
    return prepared


def best_world_match(
    reference_seq: str, world_dir: str, tie_margin: float = 1.0
) -> Tuple[Optional[str], float]:
    """
    Pick the world FASTA in ``world_dir`` whose first record best matches ``reference_seq``.

    Used to auto-map a virus to its diversity set for the conservancy step, so no manual
    virus→file table is required. Returns ``(path, identity_percent)`` or ``(None, 0.0)``.

    A single-protein reference can score ~equally against its own set and against a polyprotein
    that merely contains it (e.g. CHIKV E1 vs. the E2-E1-6K polyprotein). So among candidates
    within ``tie_margin`` percentage points of the top identity, the one whose first record is
    closest in length to the reference wins — favouring the matching single-protein set.
    """
    from .conservancy_client import max_identity  # local import: avoids import cycle at module load

    ref_len = len(reference_seq)
    candidates: List[Tuple[str, float, int]] = []  # (path, identity, first-record length)
    for fname in sorted(os.listdir(world_dir)):
        if not fname.lower().endswith((".fasta", ".fa", ".faa")):
            continue
        path = os.path.join(world_dir, fname)
        try:
            first = next(SeqIO.parse(path, "fasta"))
        except StopIteration:
            continue
        wseq = str(first.seq)
        # query the shorter into the longer so a single-protein ref maps onto a polyprotein too
        ident = max_identity(reference_seq, wseq) if ref_len <= len(wseq) else max_identity(wseq, reference_seq)
        candidates.append((path, ident, len(wseq)))

    if not candidates:
        return None, 0.0

    best_id = max(ident for _, ident, _ in candidates)
    near_top = [c for c in candidates if c[1] >= best_id - tie_margin]
    path, ident, _ = min(near_top, key=lambda c: abs(c[2] - ref_len))
    return path, ident
