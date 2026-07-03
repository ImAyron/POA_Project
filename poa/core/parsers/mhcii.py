"""Parser for IEDB MHC-II Binding Predictions (directory of HTML/TSV text outputs).

Logic preserved unchanged from the original ``mhcii.py``.
"""
from __future__ import annotations

import html
import os

import pandas as pd


def MHCIIHTMLconverter(file):
    """
    Converts a MHC-II Binding Prediction output in HTML/TSV format to a pandas DataFrame.
    """
    with open(file, "r") as raw_data:
        lines = [html.unescape(line.strip()) for line in raw_data]

    # Initialize DataFrame
    df = pd.DataFrame()
    columns = []
    for line in lines:
        # Identify header and initialize DataFrame columns
        if line.startswith("allele"):
            columns = line.split("\t")
            df = pd.DataFrame(columns=columns)
        elif line.startswith("H"):
            row_data = line.split("\t")
            # Ensure row data matches the number of columns
            row_data += ["[NaN]"] * (len(columns) - len(row_data))
            df.loc[len(df)] = row_data

    return df


def HLA_df_slice(df, HLAlist, HLAtype):
    """Filters a DataFrame for specific HLA alleles of a given class."""
    delete_alleles = [allele for allele in HLAlist if f"HLA-{HLAtype}" not in allele]
    return df[~df["allele"].isin(delete_alleles)].reset_index(drop=True)


def files_map_MHCIIBD(directory):
    """
    Scans a directory for MHCIIBP HTML result files and extracts metadata (species, protein).

    Returns:
        tuple: Lists of species, protein names, and file paths.
    """
    specie, protein, MHCII_files = [], [], []
    for _, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".html"):
                path_file = os.path.join(directory, file)
                MHCII_files.append(path_file)
                infos = file.replace(".html", "").upper().split("_")
                protein.append(infos[0])
                specie.append(infos[1])
                # Concatenate remaining parts if available for sequence ID

    return specie, protein, MHCII_files


def MHCIIAntigenEpitopes(df, specie, protein, allele_class, ic50):
    """
    Filters and formats MHC-II predicted epitopes data (NN_align by default) by IC50 and HLA type.
    """
    # Select NN_align results and add species/protein metadata
    nn_df = df[["allele", "seq_num", "start", "end", "method", "peptide",
                "smm_align_ic50", "nn_align_ic50", "nn_align_rank", "nn_align_adjusted_rank"]]
    nn_df = nn_df.assign(specie=specie, protein=protein)  # Add columns
    nn_df = nn_df[nn_df["nn_align_ic50"] != "-"].copy()
    nn_df["nn_align_ic50"] = pd.to_numeric(nn_df["nn_align_ic50"], errors="coerce")

    # Classify binding prediction by IC50 value (strong and weak bindings)
    nn_df["Prediction_NN"] = nn_df["nn_align_ic50"].apply(
        lambda x: "SB" if x < 50 else "WB" if x < 500 else "-"
    )

    # Filter for HLA allele class if specified
    if allele_class in ["DR", "DP", "DQ"]:
        nn_df = HLA_df_slice(nn_df, nn_df["allele"].unique(), allele_class)

    # Filter by IC50 threshold (NN_Align algorithm)
    nn_df = nn_df[nn_df["nn_align_ic50"] <= ic50] if ic50 > 0 else nn_df
    nn_df.reset_index(drop=True, inplace=True)

    # Construct final DataFrame
    final_df = nn_df[["specie", "protein", "allele", "start", "end", "peptide"]].copy()
    final_df.insert(0, "Method", "MHCII-Binding")
    final_df.insert(3, "ID_Sequence", "-")

    # Standardize column names
    final_df.rename(columns={
        "specie": "Specie",
        "protein": "Protein",
        "allele": "Allele",
        "start": "Initial Position",
        "end": "Final Position",
        "peptide": "Peptide Sequence",
    }, inplace=True)

    return final_df
