"""Streamlit GUI for the POA pipeline.

Launch with:  streamlit run poa/gui/app.py   (or:  python run_gui.py)

The UI logic lives in :mod:`poa.gui.app`; the UI-agnostic orchestration (working directories,
argument building, running each stage, summaries for charts) lives in :mod:`poa.gui.backend`
so it can be unit-tested without a Streamlit runtime.
"""
