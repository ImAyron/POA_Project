"""FASTA output writers.

* :func:`filesforEptConsAnalysis` — one FASTA per species for the Conservancy Analysis step
  (from the original ``prepianResultsforConservancyAnalysis.py``).
* :func:`getfastafile` — POA2 FASTA export filtered by membrane topology region
  (from the original ``POA2_v1.0.py``).

Scientific logic (headers, region selection) preserved unchanged. The only robustness change is
directory creation via ``os.makedirs(..., exist_ok=True)`` instead of the previous manual walk.
"""
from __future__ import annotations

import os


def Epitope_EptConsAnalysis(dataframe, path):
    """Creates FASTA files for Epitope Conservancy Analysis, organizing epitopes by species."""
    # List unique species in the DataFrame
    Specieslist = list(dataframe["Specie"].unique())

    # Process each species
    for sp in Specieslist:
        # Filter epitopes for the current species
        df_Specie = dataframe.loc[(dataframe["Specie"] == sp)].reset_index(drop=True)

        # Create a FASTA file for the species
        with open(f"{path}/Conservancy Analysis/{sp}_epitopes.fasta", "w") as epitopes:
            for linha in range(len(df_Specie["Peptide Sequence"])):
                # Create the FASTA header  (>Specie_Protein_Method_Init_Final  — exactly 5 fields)
                id_epitope = f">{df_Specie['Specie'][linha]}_{df_Specie['Protein'][linha]}_{df_Specie['Method'][linha]}_{df_Specie['Initial Position'][linha]}_{df_Specie['Final Position'][linha]}"
                # Write the header and sequence to the file
                epitopes.write(f"{id_epitope}\n")
                epitopes.write(f"{df_Specie['Peptide Sequence'][linha]}\n")


def filesforEptConsAnalysis(EpitopeDataframe, directory):
    """Prepares the directory and creates FASTA files for Epitope Conservancy Analysis."""
    # Ensure the "Conservancy Analysis" directory exists (robust: no failure if it already does)
    os.makedirs(f"{directory}/Conservancy Analysis", exist_ok=True)

    # Create FASTA files for epitope conservancy analysis
    Epitope_EptConsAnalysis(EpitopeDataframe, directory)


def getfastafile(path, portion_epitopes, dataframe):
    """
    Generates a FASTA file of epitopes filtered by membrane topology region.

    Parameters:
        portion_epitopes (int): 0=all, 1=outer, 2=transmembrane, 3=inner.
    """
    with open(f"{path}.fasta", "w") as results_fasta:
        for indice_line in range(len(dataframe["Epitope sequence"])):
            idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
            epitope = dataframe.iloc[indice_line]["Epitope sequence"]
            column = ""
            # Determine which membrane topology region to reference
            if portion_epitopes == 0:  # All regions
                results_fasta.write(f">{idt_epitope}")
                results_fasta.write(f"\n{epitope}\n")
            elif portion_epitopes == 1:  # Outer regions
                column = "Portion_Outside"
            elif portion_epitopes == 2:  # Transmembrane regions
                column = "Portion_TM"
            else:  # Inner regions
                column = "Portion_Inside"
            if column != "":
                if dataframe.iloc[indice_line][column] == 1:  # Belongs to the specified region
                    results_fasta.write(f">{idt_epitope}")
                    results_fasta.write(f"\n{epitope}\n")
