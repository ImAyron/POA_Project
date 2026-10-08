"""Central logging configuration for the POA pipeline.

Both the CLI and the GUI call :func:`setup_logging` once; every module gets its logger from
:func:`get_logger`. What a run prints is meant to answer three questions without reading the
code: **which stage is running**, **what it produced**, and **what failed and where**.

Every line has the same shape — ``<TAG> <stage> | <message>`` — so the log reads top to bottom
without knowing which module wrote what::

    STEP POA1/bepipred-2.0 | start
    DATA POA1/bepipred-2.0 | source: predictions.json (length filter 0-0)
    DONE POA1/bepipred-2.0 | epitopes=11 species=DENV1 elapsed=0.4s
    WARN POA2/topology     | 2 epitope(s) matched no protein in the -f FASTA
    FAIL POA2/topology     | RuntimeError: pyTMHMM is not installed elapsed=0.0s

See **message grammar** below for the tags and why they are ASCII.

Levels, and what belongs in each:

``DEBUG``
    Per-item detail — one line per file, per record, per request. Off by default; turn it on with
    ``POA_LOG_LEVEL=DEBUG`` when a stage's result needs explaining.
``INFO``
    Stage boundaries and their outcomes. One ``STEP``/``DONE`` pair per stage, with the counts
    that decide whether the run is worth continuing.
``WARNING``
    The run goes on but a result is now suspect: a service fell back to manual upload, a
    comparison set cannot discriminate, an epitope matched no protein.
``ERROR``
    The stage did not finish. Always paired with the exception, so a failure leaves a record even
    when the traceback goes somewhere else (a Streamlit page, a subprocess).

Environment:

``POA_LOG_LEVEL``
    ``DEBUG``/``INFO``/``WARNING``/``ERROR`` (default ``INFO``).
``POA_LOG_FILE``
    Path to also append the log to, at DEBUG regardless of the console level. The GUI points this
    at the analysis folder, so every run leaves ``poa.log`` next to its results.
"""
from __future__ import annotations

import logging
import os
import time
from contextlib import contextmanager
from typing import Optional

_CONFIGURED = False
_FILE_HANDLERS: dict[str, logging.Handler] = {}

#: compact for a terminal — the date is in the file log, the clock is what matters live
_CONSOLE_FORMAT = "%(asctime)s %(levelname)-7s %(name)-24s %(message)s"
_FILE_FORMAT = "%(asctime)s %(levelname)-7s %(name)-24s %(message)s"


def _level_from_env(default: int) -> int:
    name = os.environ.get("POA_LOG_LEVEL", "").strip().upper()
    return getattr(logging, name, default) if name else default


def setup_logging(level: Optional[int] = None, log_file: Optional[str] = None) -> logging.Logger:
    """
    Configure the root ``poa`` logger: one console handler, optionally one file handler.

    Idempotent for the console handler — calling it again from another entry point does not
    double every line. ``log_file`` (or ``POA_LOG_FILE``) can be added later, and is added only
    once per path, so the GUI may call this again when it binds a different analysis folder.
    """
    global _CONFIGURED
    logger = logging.getLogger("poa")
    resolved = _level_from_env(logging.INFO if level is None else level)

    if not _CONFIGURED:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(_CONSOLE_FORMAT, datefmt="%H:%M:%S"))
        handler.setLevel(resolved)
        logger.addHandler(handler)
        logger.propagate = False
        _CONFIGURED = True

    # The logger itself must pass everything the file handler may want, so it sits at the lowest
    # level in use; each handler then filters for itself.
    target = log_file or os.environ.get("POA_LOG_FILE") or None
    logger.setLevel(min(resolved, logging.DEBUG) if target else resolved)
    if target:
        use_log_file(target)
    return logger


def use_log_file(path: str) -> Optional[logging.Handler]:
    """
    Route the file log (DEBUG and up) to ``path``, **replacing** any previous file destination.

    Replacing rather than adding is the whole point: the GUI calls this every time it binds an
    analysis, and handlers that accumulated would write each line into every log file opened in
    that session — so the second analysis's lines would land in the first analysis's ``poa.log``,
    which is exactly what a per-analysis log is supposed to prevent.

    Returns the handler, or ``None`` when the file could not be opened. A log file that cannot be
    written must never stop a run — the run is the point, the record is a convenience — so the
    failure is reported on the console and execution continues with the previous destination
    left in place.
    """
    key = os.path.abspath(path)
    if key in _FILE_HANDLERS:
        return _FILE_HANDLERS[key]
    logger = logging.getLogger("poa")
    try:
        os.makedirs(os.path.dirname(key) or ".", exist_ok=True)
        handler = logging.FileHandler(key, encoding="utf-8")
    except OSError as exc:
        logger.warning("%s logging | could not open the log file %s: %s. Console only.",
                       TAG_WARN, path, exc)
        return None
    handler.setFormatter(logging.Formatter(_FILE_FORMAT, datefmt="%Y-%m-%d %H:%M:%S"))
    handler.setLevel(logging.DEBUG)

    for old_key, old in list(_FILE_HANDLERS.items()):
        logger.removeHandler(old)
        old.close()
        del _FILE_HANDLERS[old_key]

    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    _FILE_HANDLERS[key] = handler
    logger.info("%s logging | log file: %s", TAG_NOTE, key)
    return handler


def get_logger(name: str) -> logging.Logger:
    """Return a child logger under the ``poa`` namespace."""
    return logging.getLogger(f"poa.{name}")


# --------------------------------------------------------------------------- message grammar
# Every line a run produces has the same shape:
#
#     <TAG> <stage> | <message>
#
# ``TAG`` says what kind of event it is, ``stage`` says which part of the backend produced it, and
# the message carries the facts. One grammar for all of them means a log can be read top to bottom
# without knowing which module wrote what, and grepped without a regex per module:
#
#     grep 'FAIL'              -> everything that did not finish
#     grep 'POA2/'             -> the whole POA2 stage, sub-steps included
#     grep 'DONE POA1/'        -> what each POA1 method produced
#
# Tags are ASCII on purpose. The arrows and check marks used before could raise
# UnicodeEncodeError on a Windows console still running a legacy code page — the terminal this
# project is developed on — and a crash in the logging is a crash in the run.
TAG_STEP = "STEP"     # a stage started
TAG_DONE = "DONE"     # a stage finished, with what it produced
TAG_FAIL = "FAIL"     # a stage did not finish, with the exception that stopped it
TAG_NOTE = "NOTE"     # a fact established inside a stage
TAG_WARN = "WARN"     # the stage goes on, but a result is now suspect
TAG_DATA = "DATA"     # per-item detail, DEBUG only


def _format_facts(facts: dict) -> str:
    """``{'epitopes': 19}`` -> ``'epitopes=19'`` — one shape for every outcome line."""
    return " ".join(f"{key}={value}" for key, value in facts.items())


class Step:
    """
    One stage of a run, as the log sees it. Created by :func:`log_step`, not directly.

    ``stage`` is the backend path that produced the line (``POA1``, ``POA1/bepipred-2.0``,
    ``POA2/topology``, ``service/iedb-mhcii``), so the log names the part of the pipeline rather
    than the module that happens to hold the call.

    ``result()`` collects the facts for the closing line. They are reported there, not where they
    are computed, so a stage's outcome is a single line in the same shape as every other.
    """

    def __init__(self, logger: logging.Logger, stage: str):
        self._logger = logger
        self.stage = stage
        self._facts: dict[str, object] = {}
        self._started = time.perf_counter()

    @property
    def elapsed(self) -> float:
        return time.perf_counter() - self._started

    def result(self, **facts) -> "Step":
        """Record facts for the closing line, e.g. ``step.result(epitopes=19, species=2)``."""
        self._facts.update(facts)
        return self

    def detail(self, message: str, *args) -> None:
        """Per-item detail — one line per file, per record. DEBUG only."""
        self._logger.debug("%s %s | " + message, TAG_DATA, self.stage, *args)

    def note(self, message: str, *args) -> None:
        """A fact worth seeing at INFO inside the stage, without ending it."""
        self._logger.info("%s %s | " + message, TAG_NOTE, self.stage, *args)

    def warn(self, message: str, *args) -> None:
        """The stage continues, but a result is now suspect."""
        self._logger.warning("%s %s | " + message, TAG_WARN, self.stage, *args)


@contextmanager
def log_step(logger: logging.Logger, stage: str):
    """
    Bracket a backend stage with a ``STEP`` line, a ``DONE`` line, and ``FAIL`` when it raises.

    ::

        with log_step(logger, "POA1/bepipred-2.0") as step:
            ...
            step.result(epitopes=len(df), species="DENV1")

    produces::

        STEP POA1/bepipred-2.0 | start
        DONE POA1/bepipred-2.0 | epitopes=11 species=DENV1 elapsed=0.4s

    The exception is logged where the stage is known and then re-raised unchanged: callers keep
    deciding what a failure means, while the log keeps a record even when the traceback is shown
    somewhere the terminal never sees (a Streamlit page, a subprocess, a background thread).
    """
    step = Step(logger, stage)
    logger.info("%s %s | start", TAG_STEP, stage)
    try:
        yield step
    except Exception as exc:
        # The exception type is part of the message, not only of the traceback: at INFO the
        # traceback is not emitted, and 'FileNotFoundError' vs 'ValueError' is usually the whole
        # diagnosis for someone reading a log they did not produce.
        logger.error("%s %s | %s: %s elapsed=%.1fs",
                     TAG_FAIL, stage, type(exc).__name__, exc, step.elapsed)
        logger.debug("%s %s | traceback follows", TAG_DATA, stage, exc_info=True)
        raise
    facts = _format_facts(step._facts)
    logger.info("%s %s | %selapsed=%.1fs",
                TAG_DONE, stage, f"{facts} " if facts else "", step.elapsed)
