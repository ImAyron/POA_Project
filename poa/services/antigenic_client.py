"""Antigenicity prediction via the local EMBOSS ``antigenic`` program.

The original PAP/IMED web tool (imed.med.ucm.es) is offline (VPN-gated as of 2026). EMBOSS
``antigenic`` implements the identical Kolaskar & Tongaonkar (1990) method the PAP site wrapped.
Install with ``conda install -c bioconda emboss`` (on Windows: WSL2 + conda, or Docker).

We run ``antigenic ... -rformat gff`` (structured output) and convert the predicted regions to
the PAP/IMED ``.txt`` layout so the existing :mod:`poa.core.parsers.papimed` parser consumes the
result unchanged. If the binary is absent, a :class:`ServiceUnavailable` is raised so the caller
can fall back to a manually supplied ``.txt``.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from ..logging_conf import get_logger, log_step
from .base import Cache, ServiceError, ServiceResult, ServiceUnavailable, cached_call

logger = get_logger("antigenic")


def run_antigenic(fasta_path: str, minlen: int = 6, emboss_bin: str = "antigenic",
                  timeout: int = 300) -> str:
    """Run EMBOSS ``antigenic`` and return the GFF report text. Raises ServiceUnavailable if missing."""
    out_path = tempfile.NamedTemporaryFile(suffix=".gff", delete=False).name
    cmd = [emboss_bin, "-sequence", fasta_path, "-minlen", str(minlen),
           "-rformat", "gff", "-outfile", out_path, "-auto"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as exc:
        raise ServiceUnavailable(
            "EMBOSS 'antigenic' not found. Install it (conda install -c bioconda emboss) "
            "or, on Windows, run via WSL2/Docker."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise ServiceUnavailable(f"EMBOSS 'antigenic' timed out after {timeout}s.") from exc
    except subprocess.CalledProcessError as exc:
        raise ServiceError(f"EMBOSS 'antigenic' failed: {exc.stderr or exc}") from exc
    try:
        with open(out_path, encoding="utf-8") as fh:
            return fh.read()
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


def parse_gff(gff_text: str) -> Dict[str, List[Tuple[int, int, float]]]:
    """Parse an EMBOSS GFF report into {seqid: [(start, end, score), ...]}."""
    features: Dict[str, List[Tuple[int, int, float]]] = {}
    for line in gff_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        seqid = fields[0]
        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError:
            continue
        try:
            score = float(fields[5])
        except ValueError:
            score = 0.0
        features.setdefault(seqid, []).append((start, end, score))
    return features


def to_pap_txt(gff_text: str, fasta_path: str) -> str:
    """
    Convert an EMBOSS antigenic GFF report into the PAP/IMED ``.txt`` format expected by
    :func:`poa.core.parsers.papimed.PAPepitopes`, slicing peptide sequences from the proteins.
    """
    features = parse_gff(gff_text)
    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}

    blocks = []
    for seqid, seq in sequences.items():
        feats = sorted(features.get(seqid, []), key=lambda f: f[0])
        if not feats:
            continue
        lines = [f">{seqid}", "n\tStart Position\tSequence\tEnd Position"]
        for n, (start, end, _score) in enumerate(feats, start=1):
            peptide = seq[start - 1:end]
            lines.append(f"{n}\t{start}\t{peptide}\t{end}")
        blocks.append("\n".join(lines))
    return "\n".join(blocks) + ("\n" if blocks else "")


def predict(fasta_path: str, minlen: int = 6, emboss_bin: str = "antigenic",
            cache: Optional[Cache] = None, use_cache: bool = True, timeout: int = 300) -> ServiceResult:
    """Cached antigenicity prediction; result content is the raw GFF report."""
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        fasta_content = fh.read()
    params = {"minlen": minlen}
    with log_step(logger, "service/emboss-antigenic") as step:
        step.note("input=%s minlen=%s", fasta_path, minlen)
        result = cached_call(
            cache, "antigenic", "emboss-6.6", params, [fasta_content],
            lambda: run_antigenic(fasta_path, minlen, emboss_bin, timeout),
            source="local", ext="gff", use_cache=use_cache,
        )
        step.result(source=result.source, bytes=len(result.content or ""))
    return result
