"""pyTMHMM client — transmembrane topology prediction (already a local dependency).

Thin wrapper over :mod:`poa.core.topology` that reports availability and surfaces a clear
:class:`ServiceUnavailable` if pyTMHMM is not installed, so the GUI can guide the user.
"""
from __future__ import annotations

from ..core import topology
from .base import ServiceUnavailable


def is_available() -> bool:
    """Return True if pyTMHMM can be imported."""
    try:
        import pyTMHMM  # noqa: F401
        return True
    except Exception:
        return False


def predict_topology(fasta_path: str):
    """Return (ids, annotations) from pyTMHMM; raise ServiceUnavailable if not installed."""
    if not is_available():
        raise ServiceUnavailable(
            "pyTMHMM is not installed. Install it (pip install pyTMHMM) to run POA2 membrane "
            "topology analysis."
        )
    return topology.pyTMHMMpredict(fasta_path)
