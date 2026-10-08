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
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from Bio import SeqIO

from ..logging_conf import TAG_WARN, get_logger

logger = get_logger("realdata_import")

# bepipred_<virus>.json  ->  <virus>
_FNAME_RE = re.compile(r"bepipred[_-]?(?P<virus>[a-zA-Z0-9]+)", re.IGNORECASE)


@dataclass
class AntigenEntry:
    """One antigen block of a BepiPred-2.0 JSON, mapped to the pipeline header convention."""

    original_key: str
    protein: str
    specie: str
    id_seq: str
    sequence: str

    @property
    def header(self) -> str:
        return f"{self.protein}_{self.specie}_{self.id_seq}"


@dataclass
class Prepared:
    """One adapted BepiPred-2.0 JSON — one or more antigens, possibly of different species."""

    specie: str
    protein: str
    id_seq: str
    sequence: str
    reference_fasta: str  # path to the rebuilt -f FASTA
    b2_json: str          # path to the rewritten -b2 JSON
    has_x: bool           # reference contains an ambiguous 'X' (POA1 will refuse it)
    entries: List[AntigenEntry] = field(default_factory=list)

    @property
    def header(self) -> str:
        return f"{self.protein}_{self.specie}_{self.id_seq}"

    @property
    def species(self) -> List[str]:
        """Every species covered by this JSON, in file order."""
        seen: List[str] = []
        for entry in self.entries:
            if entry.specie not in seen:
                seen.append(entry.specie)
        return seen or [self.specie]


def antigen_keys(json_path: str) -> List[str]:
    """The antigen keys of a BepiPred-2.0 JSON, in file order (for building a mapping UI)."""
    with open(json_path, "r") as fh:
        data = json.load(fh)
    return list(data.get("antigens", {}))


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
    mapping: Optional[Dict[str, Tuple[str, str, str]]] = None,
) -> Prepared:
    """
    Adapt one BepiPred-2.0 JSON export to pipeline conventions.

    Writes, into ``out_dir``:
      * ``<stem>_ref.fasta`` — the reference protein(s) (``-f``), headers ``Protein_Specie_ID``.
      * ``<stem>_b2.json``   — a copy whose antigen key(s) follow ``Protein_Specie_ID`` so the
                               BepiPred-2.0 parser extracts species/protein correctly.

    Parameters:
        specie / protein / id_seq: labels applied to *every* antigen in the file.
        mapping: optional per-antigen override, ``{original_key: (protein, specie, id_seq)}``.
            Needed when a single BepiPred run covered several organisms — a submission of
            ``denv1_ns1`` + ``denv2_ns1`` comes back keyed ``denv1``/``denv2``, which the parser's
            ``Protein_Specie_ID`` regex cannot resolve, so both species would be dropped.

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
    entries: List[AntigenEntry] = []
    for i, (old_key, antigen) in enumerate(items, start=1):
        if mapping and old_key in mapping:
            prot_i, sp_i, id_i = mapping[old_key]
            id_i = id_i or id_seq
        else:
            prot_i, sp_i = protein, specie
            id_i = f"{id_seq}{i}" if multi else id_seq

        entry = AntigenEntry(old_key, prot_i, sp_i, id_i, _antigen_sequence(antigen))
        new_key, n = entry.header, 2
        while new_key in new_antigens:   # never drop an antigen to a header collision
            new_key = f"{entry.header}{n}"
            n += 1
        new_antigens[new_key] = antigen
        entries.append(entry)
    data["antigens"] = new_antigens

    species = []
    for entry in entries:
        if entry.specie not in species:
            species.append(entry.specie)
    stem = "-".join(species) if mapping else specie

    b2_path = os.path.join(out_dir, f"{stem}_b2.json")
    with open(b2_path, "w") as fh:
        json.dump(data, fh)

    ref_path = os.path.join(out_dir, f"{stem}_ref.fasta")
    with open(ref_path, "w") as fh:
        for key, entry in zip(new_antigens, entries):
            fh.write(f">{key}\n{entry.sequence}\n")

    full_seq = "".join(entry.sequence for entry in entries)
    has_x = "x" in full_seq.lower()
    if has_x:
        logger.warning("%s service/realdata-import | %s: reference for %s contains 'X' — "
                       "POA1 will refuse it (-f integrity check).", TAG_WARN, json_path, stem)

    return Prepared(
        specie=species[0] if mapping else specie,
        protein=entries[0].protein if mapping else protein,
        id_seq=id_seq,
        sequence=full_seq,
        reference_fasta=ref_path,
        b2_json=b2_path,
        has_x=has_x,
        entries=entries,
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
