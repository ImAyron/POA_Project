"""Validated paths and publication of complete batches of generated artifacts."""
from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
import uuid
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path, PureWindowsPath

HISTORY_DIRNAME = ".poa_history"

#: How many replaced batches stay on disk under ``.poa_history``.
#:
#: Nothing in the codebase ever reads an archive back — there is no restore command, no GUI
#: listing, no manifest pointing at one. The archive exists for a single manual scenario: the
#: researcher re-ran a step and goes into the results folder to fetch by hand the batch that
#: re-run replaced — including the case where it replaced it with nothing, which is how an
#: empty batch legitimately clears a directory. That scenario only ever reaches for a recent
#: generation, while every older one keeps charging full price, since a generation is a whole
#: copy of the step's output (report, reference, every spreadsheet, the entire "Conservancy
#: Analysis" tree). Keeping a few covers re-running twice before noticing; keeping all of them
#: makes the results folder grow linearly with re-runs, which only ever surfaces as a
#: surprising byte count in the delete dialog.
HISTORY_GENERATIONS = 3


def safe_child(directory, filename) -> Path:
    """Accept one filename on both Windows and POSIX, never a path or device name."""
    name = str(filename)
    if (not name or name in (".", "..") or name.rstrip(" .") != name
            or any(c in name for c in '/\\:<>"|?*')
            or any(ord(c) < 32 for c in name) or PureWindowsPath(name).is_reserved()):
        raise ValueError(f"Invalid filename: {name!r}")
    root = Path(directory).resolve()
    path = root / name
    if path.resolve().parent != root:
        raise ValueError(f"Path escapes destination: {name!r}")
    return path


def file_hash(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _prune_history(history: Path, newest: Path) -> None:
    """Drop the oldest generations, never ``newest``.

    Generation names start with a timestamp, so sorting them by name sorts them by age. The
    batch archived by the publication that just succeeded is excluded by identity rather than
    trusted to sort last: it is the only copy of what was just replaced, and the wall clock
    can step backwards (NTP, DST) between two publications.
    """
    try:
        generations = sorted(p for p in history.iterdir() if p.is_dir())
    except OSError:
        return  # no history yet, or unreadable — nothing to prune
    for stale in generations[:-HISTORY_GENERATIONS]:
        if stale != newest:
            # Cleanup is not part of the publication: the batch is already installed, and a
            # file someone still holds open must not turn a successful run into a failure.
            shutil.rmtree(stale, ignore_errors=True)


@contextmanager
def staged_outputs(directory, managed=(), publish_empty=False):
    """Publish only after all writes succeed; archive replaced outputs and roll back errors.

    Only explicitly managed names and newly generated names are replaced. Unrelated files
    stay untouched. Archives are excluded by the conservancy directory scanner.

    Args:
        directory: where the batch is published.
        managed: names this batch owns, replaced or removed even when it does not rewrite them.
        publish_empty: the body writing nothing is a result, so archive the managed outputs and
            install nothing in their place. Off by default — see below.

    A stage left empty says two things the directory cannot tell apart: "the batch is
    legitimately empty, retire the previous one" — which is how a re-run that finds no species
    clears the CSVs of species that no longer exist — and "the body gave up without writing
    anything". Read the second way, publishing archives every managed output and installs
    nothing over it: results removed because none were produced. So the caller that means the
    first one says so, and silence keeps the previous batch.
    """
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".poa_stage_", dir=root) as temp:
        stage = Path(temp)
        yield stage
        new_names = {p.name for p in stage.iterdir()}
        if not new_names and not publish_empty:
            return
        names = set(managed) | new_names
        targets = {name: safe_child(root, name) for name in names}
        if (all((stage / n).is_file() and targets[n].is_file()
                and (stage / n).read_bytes() == targets[n].read_bytes() for n in new_names)
                and not any(targets[n].exists() for n in names - new_names)):
            return
        history = root / HISTORY_DIRNAME
        # Timestamp first so the generations sort by age and a human browsing the folder can
        # tell which run each one came from; the suffix keeps two publications in the same
        # microsecond apart.
        archive = history / f"{datetime.now():%Y%m%d-%H%M%S-%f}-{uuid.uuid4().hex[:8]}"
        moved, installed = [], []
        try:
            for name in sorted(names):
                target = targets[name]
                if target.exists():
                    archive.mkdir(parents=True, exist_ok=True)
                    os.replace(target, archive / name)
                    moved.append(name)
            for name in sorted(new_names):
                os.replace(stage / name, targets[name])
                installed.append(name)
        except BaseException:
            for name in reversed(installed):
                os.replace(targets[name], stage / name)
            for name in reversed(moved):
                os.replace(archive / name, targets[name])
            raise
        _prune_history(history, archive)
