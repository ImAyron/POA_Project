"""Central logging configuration for the POA pipeline.

Both the CLI and the GUI call :func:`setup_logging` once; all modules obtain their logger
through :func:`get_logger`. This gives clear, timestamped messages and makes it obvious when
an external service failed and a manual fallback should be used.
"""
from __future__ import annotations

import logging

_CONFIGURED = False


def setup_logging(level: int = logging.INFO) -> None:
    """Configure the root ``poa`` logger with a single stream handler (idempotent)."""
    global _CONFIGURED
    if _CONFIGURED:
        return
    logger = logging.getLogger("poa")
    logger.setLevel(level)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                                            datefmt="%H:%M:%S"))
    logger.addHandler(handler)
    logger.propagate = False
    _CONFIGURED = True


def get_logger(name: str) -> logging.Logger:
    """Return a child logger under the ``poa`` namespace."""
    return logging.getLogger(f"poa.{name}")
