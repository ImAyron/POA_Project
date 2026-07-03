"""Service layer: integration with external prediction tools.

Each client tries its primary strategy (local package / REST API / subprocess / browser) and,
on failure or unavailability, raises :class:`~poa.services.base.ServiceUnavailable` so the
caller can fall back to a manually supplied result file (the original semi-automatic behavior).
Results are cached on disk (see :mod:`poa.services.cache`) to avoid re-submitting.
"""
