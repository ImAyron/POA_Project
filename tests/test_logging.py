"""The logging contract: what a run must leave behind on the terminal and on disk.

Logging is the only account of a run that survives it, and this project's workflow depends on
that: an analysis is tuned in one environment and executed in another, so the record has to travel
with the results rather than with the terminal someone happened to be looking at. These tests pin
the parts that are easy to break silently — a stage that stops reporting its counts, a failure
that leaves no line, a second ``setup_logging`` call that doubles every message.
"""
import logging
import re
from pathlib import Path

import pytest

from poa import logging_conf
from poa.logging_conf import (
    TAG_DATA, TAG_DONE, TAG_FAIL, TAG_NOTE, TAG_STEP, TAG_WARN,
    use_log_file, get_logger, log_step, setup_logging,
)


@pytest.fixture
def fresh_logging(monkeypatch):
    """A clean ``poa`` logger per test — the configuration is process-global by design."""
    logger = logging.getLogger("poa")
    saved_handlers, saved_level = list(logger.handlers), logger.level
    logger.handlers = []
    monkeypatch.setattr(logging_conf, "_CONFIGURED", False)
    monkeypatch.setattr(logging_conf, "_FILE_HANDLERS", {})
    monkeypatch.delenv("POA_LOG_LEVEL", raising=False)
    monkeypatch.delenv("POA_LOG_FILE", raising=False)
    yield logger
    logger.handlers, logger.level = saved_handlers, saved_level


# --------------------------------------------------------------------------- stage boundaries
def test_step_logs_start_and_outcome_with_its_counts(caplog):
    """A stage that only says "done" cannot be judged: 0 epitopes and 200 look the same."""
    logger = get_logger("test")
    with caplog.at_level(logging.INFO, logger="poa.test"):
        with log_step(logger, "POA1/bepipred-2.0") as step:
            step.result(epitopes=11, species="DENV1")

    messages = [r.getMessage() for r in caplog.records]
    assert messages[0] == "STEP POA1/bepipred-2.0 | start"
    assert messages[-1].startswith("DONE POA1/bepipred-2.0 | ")
    assert "epitopes=11" in messages[-1] and "species=DENV1" in messages[-1]
    assert re.search(r"elapsed=\d+\.\d+s$", messages[-1]), "o tempo deve fechar a linha"


def test_step_logs_an_error_and_reraises(caplog):
    """
    A failure must leave a line even when the traceback is shown elsewhere.

    In the GUI the exception is rendered on a Streamlit page and the terminal would otherwise show
    nothing at all; the stage is also the only place that knows *which* stage failed. The
    exception *type* is in the message because at INFO there is no traceback to read it from.
    """
    logger = get_logger("test")
    with caplog.at_level(logging.INFO, logger="poa.test"):
        with pytest.raises(RuntimeError, match="pyTMHMM"):
            with log_step(logger, "POA2/topology"):
                raise RuntimeError("pyTMHMM is not installed")

    errors = [r for r in caplog.records if r.levelno == logging.ERROR]
    assert len(errors) == 1
    message = errors[0].getMessage()
    assert message.startswith("FAIL POA2/topology | ")
    assert "RuntimeError: pyTMHMM is not installed" in message
    assert re.search(r"elapsed=\d+\.\d+s$", message)


def test_step_detail_is_debug_only(caplog):
    """Per-item lines must not drown the stage boundaries at the default level."""
    logger = get_logger("test")
    with caplog.at_level(logging.INFO, logger="poa.test"):
        with log_step(logger, "stage") as step:
            step.detail("one line per file")
            step.note("worth seeing")

    messages = [r.getMessage() for r in caplog.records]
    assert not any("one line per file" in m for m in messages)
    assert "NOTE stage | worth seeing" in messages


def test_step_without_facts_still_closes_cleanly(caplog):
    """Not every stage has a count; it must still close rather than print 'None'."""
    logger = get_logger("test")
    with caplog.at_level(logging.INFO, logger="poa.test"):
        with log_step(logger, "stage"):
            pass

    closing = [r.getMessage() for r in caplog.records][-1]
    assert "None" not in closing
    assert re.fullmatch(r"DONE stage \| elapsed=\d+\.\d+s", closing)


def test_every_line_is_tagged_and_names_its_stage(caplog):
    """
    One grammar for every line: ``<TAG> <stage> | <message>``.

    This is what lets a log be read without knowing which module wrote what, and grepped by stage
    (``POA2/``) or by outcome (``FAIL``) instead of by a regex per module.
    """
    logger = get_logger("test")
    with caplog.at_level(logging.DEBUG, logger="poa.test"):
        # a stage that succeeds covers STEP/NOTE/DATA/WARN/DONE
        with log_step(logger, "POA2/conservancy") as step:
            step.note("criterion=%s", ">=70")
            step.detail("row %d", 1)
            step.warn("nothing survived")
        # and one that raises covers FAIL — DONE is never reached on that path
        try:
            with log_step(logger, "POA2/conservancy"):
                raise ValueError("boom")
        except ValueError:
            pass

    messages = [r.getMessage() for r in caplog.records]
    for message in messages:
        assert re.match(r"^\w+ POA2/conservancy \| .", message), message
    tags = {m.split(" ", 1)[0] for m in messages}
    assert tags == {TAG_STEP, TAG_DONE, TAG_FAIL, TAG_NOTE, TAG_WARN, TAG_DATA}


def test_tags_are_ascii_only():
    """
    A Windows console on a legacy code page raises UnicodeEncodeError on non-ASCII output.

    This project is developed on Windows, and a crash inside the logging is a crash in the run —
    which is why the arrows and check marks used at first were replaced by these tags.
    """
    for tag in (TAG_STEP, TAG_DONE, TAG_FAIL, TAG_NOTE, TAG_WARN, TAG_DATA):
        tag.encode("ascii")           # raises if it ever stops being ASCII


# --------------------------------------------------------------------------- configuration
def test_setup_logging_is_idempotent(fresh_logging):
    """Both the CLI and the GUI call it; a second call used to be harmless and must stay so."""
    setup_logging()
    setup_logging()
    streams = [h for h in fresh_logging.handlers if isinstance(h, logging.StreamHandler)]
    assert len(streams) == 1, "uma segunda chamada duplicaria cada linha da saída"


def test_log_level_comes_from_the_environment(fresh_logging, monkeypatch):
    monkeypatch.setenv("POA_LOG_LEVEL", "DEBUG")
    setup_logging()
    assert fresh_logging.level == logging.DEBUG


def test_invalid_log_level_falls_back_instead_of_crashing(fresh_logging, monkeypatch):
    """A typo in an env var must not stop a run."""
    monkeypatch.setenv("POA_LOG_LEVEL", "VERBOSE")
    setup_logging()
    assert fresh_logging.level == logging.INFO


# --------------------------------------------------------------------------- the file log
def test_log_file_records_debug_even_when_the_console_is_at_info(fresh_logging, tmp_path):
    """
    The file is the record; the console is the live view.

    The console stays readable at INFO while the file keeps the per-item detail, so a result can
    be explained after the fact without re-running anything.
    """
    target = tmp_path / "analysis" / "poa.log"
    setup_logging(log_file=str(target))
    logger = get_logger("test")
    with log_step(logger, "stage") as step:
        step.detail("per-item detail %d", 42)

    written = target.read_text(encoding="utf-8")
    assert "STEP stage | start" in written
    assert "DATA stage | per-item detail 42" in written      # DEBUG reached the file
    assert "DONE stage | elapsed=" in written


def test_log_file_is_added_once_per_path(fresh_logging, tmp_path):
    """The GUI calls use_log_file on every rerun; the file must not get N copies of each line."""
    target = tmp_path / "poa.log"
    use_log_file(str(target))
    use_log_file(str(target))
    use_log_file(str(tmp_path / "." / "poa.log"))     # same file, different spelling

    files = [h for h in fresh_logging.handlers if isinstance(h, logging.FileHandler)]
    assert len(files) == 1


def test_switching_the_log_file_stops_writing_to_the_previous_one(fresh_logging, tmp_path):
    """
    One analysis's lines must never land in another analysis's log.

    The GUI calls this whenever it binds an analysis. While the handlers accumulated, every line
    was written to every file opened in the session, so opening a second analysis appended its
    whole run to the first one's ``poa.log`` — silently undoing the one thing a per-analysis log
    is for.
    """
    first, second = tmp_path / "A" / "poa.log", tmp_path / "B" / "poa.log"
    setup_logging()
    logger = get_logger("test")

    use_log_file(str(first))
    logger.info("belongs to A")
    use_log_file(str(second))
    logger.info("belongs to B")

    assert "belongs to A" in first.read_text(encoding="utf-8")
    assert "belongs to B" not in first.read_text(encoding="utf-8")
    assert "belongs to B" in second.read_text(encoding="utf-8")
    files = [h for h in fresh_logging.handlers if isinstance(h, logging.FileHandler)]
    assert len(files) == 1, "o handler anterior precisa ser removido, não apenas somado"


def test_unwritable_log_file_warns_and_lets_the_run_continue(fresh_logging, tmp_path, caplog):
    """The run is the point and the record is a convenience — never the other way round."""
    setup_logging()
    blocker = tmp_path / "not_a_dir"
    blocker.write_text("", encoding="utf-8")

    # setup_logging stops propagation so the app owns its output; caplog reads the root handler,
    # so it has to be let through again just for the assertion.
    fresh_logging.propagate = True
    with caplog.at_level(logging.WARNING, logger="poa"):
        handler = use_log_file(str(blocker / "poa.log"))

    assert handler is None
    assert any("log file" in r.getMessage() for r in caplog.records)


# --------------------------------------------------------------------------- structure
def test_no_two_modules_share_a_logger_name():
    """
    Two modules under one name make their lines indistinguishable.

    ``poa/core/parsers/conservancy.py`` and ``poa/services/conservancy_client.py`` both used
    ``get_logger("conservancy")``, so a message about filtering CSVs and one about computing them
    were reported by the same logger.
    """
    root = Path(__file__).resolve().parents[1] / "poa"
    seen = {}
    for path in sorted(root.rglob("*.py")):
        for name in re.findall(r'get_logger\("(\w+)"\)', path.read_text(encoding="utf-8")):
            seen.setdefault(name, []).append(path.name)

    shared = {name: files for name, files in seen.items() if len(files) > 1}
    assert not shared, f"nomes de logger repetidos: {shared}"


def test_pipeline_stages_are_bracketed_by_log_step():
    """
    Every POA stage reports its own boundary.

    The stage names are what make the log readable per step, which is the whole point of the
    instrumentation; a stage added later without a bracket would log nothing of its own.
    """
    source = (Path(__file__).resolve().parents[1] / "poa" / "core" / "pipeline.py").read_text(
        encoding="utf-8")
    stages = set(re.findall(r'log_step\(logger,\s*"([^"]+)"', source))

    assert {"POA1", "POA2", "POA2/conservancy", "POA2/topology"} <= stages
    for method in ("bepipred-2.0", "bepipred-3.0", "pap-imed", "netctl-1.2", "mhc-ii",
                   "other-predictors"):
        assert f"POA1/{method}" in stages, f"etapa sem log_step: {method}"


