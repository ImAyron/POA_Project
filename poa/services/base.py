"""Common types and helpers for the service layer: results, errors, cache glue, manual fallback."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence

from .cache import Cache, cache_key


class ServiceError(Exception):
    """A service failed in a way the caller cannot recover from automatically."""


class ServiceUnavailable(ServiceError):
    """The service (network/tool/package) is unavailable — the caller should offer manual upload."""


@dataclass
class ServiceResult:
    """Outcome of a prediction request."""

    content: str            # raw result text (TSV / CSV / FASTA / report)
    source: str             # 'cache' | 'api' | 'local' | 'browser' | 'manual'
    path: Optional[str] = None


def read_manual(file_path: str | Path) -> ServiceResult:
    """Fallback path: use a result file the user obtained/uploaded manually."""
    p = Path(file_path)
    if not p.exists():
        raise ServiceError(f"Manual result file not found: {p}")
    return ServiceResult(content=p.read_text(encoding="utf-8", errors="replace"),
                         source="manual", path=str(p))


def cached_call(
    cache: Optional[Cache],
    tool: str,
    version,
    params: dict,
    inputs: Sequence[str],
    fn: Callable[[], str],
    *,
    source: str = "live",
    ext: str = "dat",
    use_cache: bool = True,
) -> ServiceResult:
    """
    Return a cached result if present, otherwise call ``fn`` (which performs the real request)
    and cache its output. ``fn`` may raise :class:`ServiceUnavailable` / :class:`ServiceError`.
    """
    if cache is not None and use_cache:
        key = cache_key(tool, version, params, *inputs)
        cached = cache.get(tool, key, ext)
        if cached is not None:
            return ServiceResult(content=cached, source="cache", path=str(cache.path(tool, key, ext)))

    content = fn()

    if cache is not None:
        key = cache_key(tool, version, params, *inputs)
        path = cache.put(tool, key, content, ext)
        return ServiceResult(content=content, source=source, path=str(path))
    return ServiceResult(content=content, source=source)
