"""Tests for the browser-automation client (availability / graceful fallback only).

The live DTU submission is not exercised here (it depends on external services and installed
browser binaries); only the graceful-unavailability contract is verified.
"""
import pytest

from poa.services import browser_client
from poa.services.base import ServiceUnavailable


def test_is_available_returns_bool():
    assert isinstance(browser_client.is_available(), bool)


def test_netctl_without_playwright_falls_back(monkeypatch):
    # Force the "Playwright missing" path regardless of the environment.
    def boom():
        raise ServiceUnavailable("no playwright")

    monkeypatch.setattr(browser_client, "_require_playwright", boom)
    with pytest.raises(ServiceUnavailable):
        browser_client.run_netctl(">SARS_SPIKE_NP1\nMKTAYIAMK\n")


def test_bepipred2_without_playwright_falls_back(monkeypatch):
    def boom():
        raise ServiceUnavailable("no playwright")

    monkeypatch.setattr(browser_client, "_require_playwright", boom)
    with pytest.raises(ServiceUnavailable):
        browser_client.run_bepipred2(">SARS_SPIKE_NP1\nMKTAYIAMK\n")
