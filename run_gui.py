#!/usr/bin/env python
"""Launch the POA Streamlit GUI.

Equivalent to:  streamlit run poa/gui/app.py
Any extra CLI args are forwarded to Streamlit (e.g. --server.port 8502).
"""
import subprocess
import sys
from pathlib import Path


def main():
    app = Path(__file__).parent / "poa" / "gui" / "app.py"
    cmd = [sys.executable, "-m", "streamlit", "run", str(app), *sys.argv[1:]]
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
