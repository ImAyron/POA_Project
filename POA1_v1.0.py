#!/usr/bin/env python
"""POA1 — Pipeline de Otimização de Antígenos, stage 1.

Backward-compatible entry point. The implementation now lives in the ``poa`` package
(``poa.cli.poa1`` / ``poa.core.pipeline``); command-line usage is unchanged, e.g.:

    python POA1_v1.0.py -f proteins.fasta -d results -b3 bepipred3.fasta -e y
"""
from poa.cli.poa1 import main

if __name__ == "__main__":
    main()
