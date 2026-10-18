"""Parser for NetCTL 1.2 cytotoxic T-cell epitope predictions (HTML output).

Logic preserved unchanged from the original ``netctl.py``.
"""
from __future__ import annotations

import html

import pandas as pd


def netctlAntigenEpitopes(file):
    """
    Extracts predicted epitopes from NetCTL-1.2 Server results in HTML format,
    creating a DataFrame with standardized columns.
    """
    # Define column names according to NetCTL result table and initialize DataFrame
    columns = [
        "Residue_number", "ID", "Protein_identifier", "pep", "Peptide_sequence",
        "aff", "Predicted_MHC_binding_affinity", "aff_rescale", "Rescale_binding_affinity",
        "cle", "C_terminal_cleavage_affinity", "tap", "TAP_transport_efficiency",
        "COMB", "Prediction_score", "Identified_MHC_ligands",
    ]
    netctl_results_df = pd.DataFrame(columns=columns)

    # Read and parse HTML file contents
    with open(file, "r") as raw_data:
        netctl_results = [html.unescape(line.strip()) for line in raw_data if not line.isspace()]

    # Fill DataFrame with parsed NetCTL results.
    # FIX (user-approved): use a sequential row index — the file line number is non-sequential and
    # pandas>=2 rejects `.loc[i] = list` for it — and pad only SHORT rows up to the column count.
    # Epitope lines carry the trailing '<-E' marker (16 tokens); non-epitope lines lack it
    # (15 tokens) and get a '-' placeholder in the Identified_MHC_ligands column.
    row_index = 0
    for line in netctl_results:
        if line and line[0].isdigit():
            line = " ".join(line.split())  # Remove extra spaces
            new_df_line = line.split(" ")

            while len(new_df_line) < len(columns):
                new_df_line.append("-")
            if len(new_df_line) != len(columns):
                continue  # skip unexpectedly long / malformed lines

            netctl_results_df.loc[row_index] = new_df_line
            row_index += 1

    # Select predicted epitopes identified by NetCTL
    netctl_epitopes = netctl_results_df[netctl_results_df["Identified_MHC_ligands"] == "<-E"].reset_index(drop=True)

    result_columns = ["Method", "Specie", "Protein", "ID_Sequence", "Initial Position", "Final Position", "Peptide Sequence"]
    if netctl_epitopes.empty:
        return pd.DataFrame(columns=result_columns)

    # Select and rename relevant columns
    netctl_epitopes_slice = netctl_epitopes[["Protein_identifier", "Residue_number", "Peptide_sequence"]].copy()
    netctl_epitopes_slice.rename(columns={"Residue_number": "Initial Position", "Peptide_sequence": "Peptide Sequence"}, inplace=True)

    # Add standard columns
    netctl_epitopes_slice["Method"] = "NetCTL"
    netctl_epitopes_slice["ID_Sequence"] = "-"

    # Split 'Protein_identifier' to extract 'Specie' and 'Protein'
    netctl_epitopes_slice["Specie"] = netctl_epitopes_slice["Protein_identifier"].apply(lambda x: x.split("_")[1])
    netctl_epitopes_slice["Protein"] = netctl_epitopes_slice["Protein_identifier"].apply(lambda x: x.split("_")[0])

    # Calculate 'Final Position' based on 'Initial Position' and peptide length
    netctl_epitopes_slice["Final Position"] = netctl_epitopes_slice.apply(
        lambda row: int(row["Initial Position"]) + len(row["Peptide Sequence"]) - 1, axis=1
    )

    # Organize columns in final DataFrame
    results_df = netctl_epitopes_slice[["Method", "Specie", "Protein", "ID_Sequence", "Initial Position", "Final Position", "Peptide Sequence"]]

    return results_df
