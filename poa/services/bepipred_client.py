"""BepiPred client — local B-cell epitope prediction.

BepiPred-3.0 ships as an installable package (``pip install bp3`` / the ``UberClifford/BepiPred-3.0``
repo) that produces a FASTA where uppercase residues mark epitopes — directly compatible with
:func:`poa.core.parsers.bepipred.bp3_FastaAnalysis`. It is heavy (PyTorch + ESM-2), so it is an
optional, on-demand dependency.

BepiPred-2.0 has no public API/PyPI package (DTU academic download only) and its JSON output is
therefore handled via the manual-upload fallback; there is no automated 2.0 client here.

If the BepiPred-3.0 CLI is not available, a :class:`ServiceUnavailable` is raised so the caller
can fall back to a manually downloaded result file.
"""
from __future__ import annotations

import glob
import os
import subprocess
import tempfile
from typing import List, Optional

from ..logging_conf import get_logger, log_step
from .base import Cache, ServiceError, ServiceResult, ServiceUnavailable, cached_call

logger = get_logger("bepipred")

def _resolve_cli(explicit: Optional[List[str]] = None) -> List[str]:
    """Locate a runnable BepiPred-3.0 CLI, or raise ServiceUnavailable."""
    if explicit:
        return explicit
    import importlib.util
    import shutil

    if shutil.which("bepipred3_CLI.py"):
        return ["bepipred3_CLI.py"]
    if shutil.which("bp3"):
        return ["bp3"]
    if importlib.util.find_spec("bp3") is not None:
        return ["python", "-m", "bp3"]
    raise ServiceUnavailable(
        "BepiPred-3.0 CLI not found. Install it (pip install bp3) or clone "
        "github.com/UberClifford/BepiPred-3.0; otherwise upload the result FASTA manually."
    )


def run_bepipred3(fasta_path: str, out_dir: Optional[str] = None, pred: str = "vt_pred",
                  cli: Optional[List[str]] = None, timeout: int = 3600) -> str:
    """
    Run the BepiPred-3.0 CLI and return the content of its output FASTA (uppercase = epitope).
    Raises ServiceUnavailable if the CLI is missing.
    """
    command = _resolve_cli(cli)
    out_dir = out_dir or tempfile.mkdtemp(prefix="bp3_")
    os.makedirs(out_dir, exist_ok=True)
    cmd = command + ["-i", fasta_path, "-o", out_dir, "-pred", pred]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as exc:
        raise ServiceUnavailable(f"BepiPred-3.0 CLI not runnable: {exc}") from exc
    except subprocess.TimeoutExpired as exc:
        raise ServiceUnavailable(f"BepiPred-3.0 timed out after {timeout}s.") from exc
    except subprocess.CalledProcessError as exc:
        raise ServiceError(f"BepiPred-3.0 failed: {exc.stderr or exc}") from exc

    fastas = sorted(glob.glob(os.path.join(out_dir, "*.fasta")))
    if not fastas:
        raise ServiceError(f"BepiPred-3.0 produced no FASTA output in {out_dir}.")
    with open(fastas[0], encoding="utf-8") as fh:
        return fh.read()


def predict(fasta_path: str, pred: str = "vt_pred", cli: Optional[List[str]] = None,
            cache: Optional[Cache] = None, use_cache: bool = True, timeout: int = 3600) -> ServiceResult:
    """Cached BepiPred-3.0 prediction; content is the output FASTA (uppercase = epitope)."""
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        fasta_content = fh.read()
    params = {"pred": pred}
    with log_step(logger, "service/bepipred-3.0") as step:
        step.note("input=%s pred=%s", fasta_path, pred)
        result = cached_call(
            cache, "bepipred3", "bp3", params, [fasta_content],
            lambda: run_bepipred3(fasta_path, None, pred, cli, timeout),
            source="local", ext="fasta", use_cache=use_cache,
        )
        step.result(source=result.source, bytes=len(result.content or ""))
    return result
