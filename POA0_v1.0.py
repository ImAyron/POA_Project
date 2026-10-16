#!/usr/bin/env python
"""POA0 — Pipeline de Otimização de Antígenos, stage 0: database search.

Entry point in the same style as POA1/POA2. The implementation lives in the ``poa`` package
(``poa.cli.poa0`` / ``poa.core.dbsearch``), e.g.:

    python POA0_v1.0.py -s "Oropouche virus" -l proteinas
    python POA0_v1.0.py -a A0A0D4BSW3 -r 482-703 --nome GcHead --sigla OROV -d entradas
"""
import os
import sys

# With PYTHONSAFEPATH=1 the script's own directory is not added to sys.path, so the sibling
# ``poa`` package would not import.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from poa.cli.poa0 import main

if __name__ == "__main__":
    sys.exit(main())
