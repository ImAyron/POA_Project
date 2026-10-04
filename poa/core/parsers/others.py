"""Parser for arbitrary predictors supplied as a standardized FASTA file.

Preserves the original length filters and supports wrapped FASTA records.
"""
from __future__ import annotations

import pandas as pd
from Bio import SeqIO


def fasta_epitopes(fasta, lenght_min, lenght_max):
    """
    Extracts and organizes epitopes from a FASTA file, applying optional length constraints.

    Header format:
        >Protein_Specie_Method_[NP_]ID_Init_Final
    """
    # Create a new DataFrame to store results
    results_df = pd.DataFrame(columns=["Method", "Specie", "Protein", "ID_Sequence", "Initial Position", "Final Position", "Peptide Sequence"])

    # A FASTA record may span several lines; line wrapping must not split an epitope.
    for index, record in enumerate(SeqIO.parse(fasta, "fasta")):
        fields = record.id.upper().split("_")
        if len(fields) < 6 or not all(fields):
            raise ValueError(
                f"Invalid epitope header '{record.id}'; expected Protein_Specie_Method_ID_Init_Final."
            )
        protein, species, method = fields[:3]
        results_df.loc[index] = [
            method, species, protein, "_".join(fields[3:-2]),
            fields[-2], fields[-1], str(record.seq).upper(),
        ]

    # Apply length constraints to epitopes
    if lenght_min == 0 and lenght_max == 0:
        pass  # No length filtering
    elif lenght_min == 0:
        results_df = results_df[results_df["Peptide Sequence"].str.len() <= lenght_max]
    elif lenght_max == 0:
        results_df = results_df[results_df["Peptide Sequence"].str.len() >= lenght_min]
    else:
        results_df = results_df[
            (results_df["Peptide Sequence"].str.len() >= lenght_min) &
            (results_df["Peptide Sequence"].str.len() <= lenght_max)
        ]

    results_df = results_df.reset_index(drop=True)

    return results_df
