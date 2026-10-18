"""On-disk cache for external prediction results, keyed by a hash of tool + params + input.

Avoids re-submitting the same sequences to the same tool/parameters. Cache entries are plain
files under ``base_dir/<tool>/<key>.<ext>`` so they are easy to inspect and to clear.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Optional

DEFAULT_CACHE_DIR = ".poa_cache"


def cache_key(tool: str, version, params: dict, *input_contents) -> str:
    """Deterministic SHA-256 key from the tool, its version, parameters and input content(s)."""
    h = hashlib.sha256()

    def _feed(value):
        if isinstance(value, bytes):
            h.update(value)
        else:
            h.update(str(value).encode("utf-8"))
        h.update(b"\x00")

    _feed(tool)
    _feed(version)
    _feed(json.dumps(params, sort_keys=True, default=str))
    for content in input_contents:
        _feed(content)
    return h.hexdigest()


class Cache:
    """Simple file-backed cache. Pass ``None`` where a cache is optional to disable it."""

    def __init__(self, base_dir: str | Path = DEFAULT_CACHE_DIR):
        self.base = Path(base_dir)

    def _path(self, tool: str, key: str, ext: str) -> Path:
        return self.base / tool / f"{key}.{ext}"

    def path(self, tool: str, key: str, ext: str = "dat") -> Path:
        return self._path(tool, key, ext)

    def get(self, tool: str, key: str, ext: str = "dat") -> Optional[str]:
        p = self._path(tool, key, ext)
        if p.exists():
            return p.read_text(encoding="utf-8")
        return None

    def put(self, tool: str, key: str, content: str, ext: str = "dat") -> Path:
        p = self._path(tool, key, ext)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content, encoding="utf-8")
        return p

    def clear(self, tool: Optional[str] = None) -> None:
        import shutil

        target = self.base / tool if tool else self.base
        if target.exists():
            shutil.rmtree(target)
