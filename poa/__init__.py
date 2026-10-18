"""POA — Pipeline de Otimização de Antígenos (Antigen Optimization Pipeline).

Refactored package layout:

* ``poa.core``     — parsing + pipeline logic (pure, no network). Reusable by CLI and GUI.
* ``poa.services`` — integration with external prediction services (network / subprocess).
* ``poa.cli``      — command-line entry points (backward compatible with the original scripts).
* ``poa.gui``      — Streamlit interface.

The original scientific ranking/selection logic is preserved unchanged; this package only
reorganizes it into layers and adds error handling, logging, caching and tests around it.
"""

__version__ = "1.0.0"
