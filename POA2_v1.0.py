#!/usr/bin/env python
"""POA2 — Pipeline de Otimização de Antígenos, stage 2.

Backward-compatible entry point. The implementation now lives in the ``poa`` package
(``poa.cli.poa2`` / ``poa.core.pipeline``); command-line usage is unchanged, e.g.:

    python POA2_v1.0.py -g True -t 70 -d conservancy_csvs -f proteins.fasta -rf 1
"""
from poa.cli.poa2 import main

if __name__ == "__main__":
    main()
