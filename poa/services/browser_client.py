"""Best-effort browser automation for tools with no API and no Windows binary.

Covers NetCTL 1.2 and BepiPred-2.0 (DTU Health Tech), which expose only asynchronous web forms
(submit -> job id -> poll -> results). This is inherently fragile: the DTU pages and their result
markup can change, so treat this as *best effort* with the manual upload always available as the
reliable fallback (the GUI defaults these tools to manual).

Requires Playwright:
    pip install playwright
    playwright install chromium

EXPERIMENTAL: the form/result selectors below are conservative defaults and may need tuning
against the live pages. Every failure path raises :class:`ServiceUnavailable`/:class:`ServiceError`
so the caller cleanly offers manual upload instead.
"""
from __future__ import annotations

import re
from typing import Callable, Optional

from ..logging_conf import TAG_NOTE, get_logger
from .base import Cache, ServiceError, ServiceResult, ServiceUnavailable, cached_call

logger = get_logger("browser")

DTU_NETCTL_URL = "https://services.healthtech.dtu.dk/services/NetCTL-1.2/"
DTU_BEPIPRED2_URL = "https://services.healthtech.dtu.dk/services/BepiPred-2.0/"


def is_available() -> bool:
    """Return True if Playwright (Python API) is importable."""
    try:
        import playwright.sync_api  # noqa: F401
        return True
    except Exception:
        return False


def _require_playwright():
    try:
        from playwright.sync_api import sync_playwright
        return sync_playwright
    except ImportError as exc:
        raise ServiceUnavailable(
            "Playwright is not installed. Run: pip install playwright && playwright install chromium. "
            "Meanwhile, submit the sequences on the tool's website and upload the result file (manual fallback)."
        ) from exc


def submit_dtu_form(
    url: str,
    fasta_text: str,
    *,
    extra_fill: Optional[Callable] = None,
    timeout_ms: int = 240000,
    headless: bool = True,
) -> str:
    """
    Submit a FASTA to a DTU Health Tech service form and return the results page text.

    Parameters:
        url: the service URL (e.g. :data:`DTU_NETCTL_URL`).
        fasta_text: FASTA content to paste into the sequence textarea.
        extra_fill: optional callable(page) to set service-specific options before submitting.
        timeout_ms: overall navigation/wait timeout.
        headless: run the browser headless.

    Raises:
        ServiceUnavailable: Playwright/browser missing or the page could not be driven.
        ServiceError: submission ran but no result could be extracted.
    """
    sync_playwright = _require_playwright()
    logger.info("%s service/browser | submitting to %s", TAG_NOTE, url)
    try:
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch(headless=headless)
            except Exception as exc:  # browser binaries not installed
                raise ServiceUnavailable(
                    f"Could not launch Chromium ({exc}). Run: playwright install chromium."
                ) from exc
            try:
                page = browser.new_page()
                page.set_default_timeout(timeout_ms)
                page.goto(url, wait_until="domcontentloaded")

                textarea = page.locator("textarea").first
                textarea.wait_for(state="visible", timeout=timeout_ms)
                textarea.fill(fasta_text)

                if extra_fill is not None:
                    extra_fill(page)

                # Best-effort submit: a submit button/input.
                submitted = False
                for getter in (
                    lambda: page.get_by_role("button", name=re.compile("submit", re.I)).first,
                    lambda: page.locator("input[type=submit]").first,
                ):
                    try:
                        getter().click()
                        submitted = True
                        break
                    except Exception:
                        continue
                if not submitted:
                    raise ServiceError("Could not locate a submit control on the DTU form.")

                # DTU shows a "job in progress" page that refreshes to the results.
                page.wait_for_load_state("networkidle", timeout=timeout_ms)
                content = page.content()
                return content
            finally:
                browser.close()
    except (ServiceUnavailable, ServiceError):
        raise
    except Exception as exc:  # any other Playwright/runtime failure -> unavailable, fall back to manual
        raise ServiceUnavailable(f"Browser automation failed ({exc}). Use manual upload as fallback.") from exc


def run_netctl(fasta_text: str, cache: Optional[Cache] = None, use_cache: bool = True,
               headless: bool = True) -> ServiceResult:
    """Best-effort NetCTL 1.2 submission. Returns the raw results HTML (feed to the NetCTL parser)."""
    return cached_call(
        cache, "netctl", "dtu-web", {}, [fasta_text],
        lambda: submit_dtu_form(DTU_NETCTL_URL, fasta_text, headless=headless),
        source="browser", ext="html", use_cache=use_cache,
    )


def run_bepipred2(fasta_text: str, cache: Optional[Cache] = None, use_cache: bool = True,
                  headless: bool = True) -> ServiceResult:
    """Best-effort BepiPred-2.0 submission. Returns the raw results HTML/JSON page text."""
    return cached_call(
        cache, "bepipred2", "dtu-web", {}, [fasta_text],
        lambda: submit_dtu_form(DTU_BEPIPRED2_URL, fasta_text, headless=headless),
        source="browser", ext="html", use_cache=use_cache,
    )
