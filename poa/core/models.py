"""Shared data-model constants for the POA pipeline.

Every POA1 parser produces a DataFrame with the same standardized column set, so the
consolidation step in :mod:`poa.core.pipeline` can concatenate them uniformly.
"""
from __future__ import annotations

#: Standardized columns produced by every POA1 epitope parser.
STANDARD_COLUMNS = [
    "Method",
    "Specie",
    "Protein",
    "ID_Sequence",
    "Initial Position",
    "Final Position",
    "Peptide Sequence",
]

#: The MHC-II parser additionally carries the HLA allele before consolidation.
MHCII_COLUMNS = [
    "Method",
    "Specie",
    "Protein",
    "Allele",
    "ID_Sequence",
    "Initial Position",
    "Final Position",
    "Peptide Sequence",
]
