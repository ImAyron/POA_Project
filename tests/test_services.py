"""Tests for the service-layer scaffolding: cache and manual fallback."""
import pytest

from poa.services import base
from poa.services.cache import Cache, cache_key


def test_cache_key_is_deterministic_and_sensitive():
    k1 = cache_key("mhcii", "1", {"allele": "DR"}, "SEQ")
    k2 = cache_key("mhcii", "1", {"allele": "DR"}, "SEQ")
    k3 = cache_key("mhcii", "1", {"allele": "DP"}, "SEQ")
    k4 = cache_key("mhcii", "1", {"allele": "DR"}, "OTHER")
    assert k1 == k2
    assert k1 != k3  # different params
    assert k1 != k4  # different input


def test_cache_put_get_roundtrip(tmp_path):
    cache = Cache(tmp_path)
    assert cache.get("tool", "abc") is None
    cache.put("tool", "abc", "hello")
    assert cache.get("tool", "abc") == "hello"


def test_read_manual(tmp_path):
    f = tmp_path / "result.tsv"
    f.write_text("allele\tpeptide\n")
    res = base.read_manual(str(f))
    assert res.source == "manual"
    assert "allele" in res.content


def test_read_manual_missing_raises(tmp_path):
    with pytest.raises(base.ServiceError):
        base.read_manual(str(tmp_path / "nope.tsv"))


def test_cached_call_caches_and_calls_once(tmp_path):
    cache = Cache(tmp_path)
    calls = {"n": 0}

    def fn():
        calls["n"] += 1
        return "RESULT"

    r1 = base.cached_call(cache, "tool", "1", {"p": 1}, ["input"], fn, source="api")
    assert r1.source == "api" and r1.content == "RESULT"

    r2 = base.cached_call(cache, "tool", "1", {"p": 1}, ["input"], fn, source="api")
    assert r2.source == "cache" and r2.content == "RESULT"
    assert calls["n"] == 1  # fn not called again


def test_cached_call_propagates_unavailable(tmp_path):
    cache = Cache(tmp_path)

    def fn():
        raise base.ServiceUnavailable("network down")

    with pytest.raises(base.ServiceUnavailable):
        base.cached_call(cache, "tool", "1", {}, ["x"], fn)
