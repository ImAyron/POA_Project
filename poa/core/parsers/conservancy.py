"""Parser and filter for IEDB Epitope Conservancy Analysis results (CSV files).

Logic preserved unchanged from the original ``conservancyAnalysis.py``. This module reads and
filters the CSVs; the *local reimplementation* that generates such results lives in
:mod:`poa.services.conservancy_client` (added later).
"""
from __future__ import annotations

import os

import pandas as pd


def map_EpConservFiles(directory):
    """Retrieve file paths for all CSV files in ``directory`` and subdirectories."""
    EpConservFiles = []
    for folders, subfolders, files in os.walk(directory):
        for file in files:
            path_file = os.path.join(folders, file)
            path_file = str(path_file)
            if ".csv" in file:
                EpConservFiles.append(path_file)
    return EpConservFiles


def EpitConservAnalysis(ID_threshold, symbol, seq_match, max_ID, min_ID, directory):
    """
    Processes and filters the results of Epitope Conservancy Analysis from multiple CSV files.

    Returns:
        pd.DataFrame: filtered and processed conservancy analysis results.
    """
    # Get the list of CSV files in the directory
    filesList = map_EpConservFiles(directory)

    # Initialize an empty DataFrame to store consolidated results
    allSp_ConservDF = pd.DataFrame(columns=["Epitope #", "Epitope name", "Epitope sequence", "Epitope length", "Percent of protein sequence matches at identity <= 100%", "Minimum identity", "Maximum identity", "View details"])

    # Consolidate data from all CSV files
    for file in filesList:
        Sp_ConservDF = pd.read_csv(file)  # Read each CSV file
        allSp_ConservDF = pd.merge(allSp_ConservDF, Sp_ConservDF, how="outer")  # Merge with the main dataframe

    # Clean and preprocess the data
    allSp_ConservDF["Percent of protein sequence matches at identity <= 100%"] = allSp_ConservDF["Percent of protein sequence matches at identity <= 100%"].map(lambda x: x.replace("%", ""))
    allSp_ConservDF["Minimum identity"] = allSp_ConservDF["Minimum identity"].map(lambda x: x.lstrip("").rstrip("%"))  # removing the '%' character
    allSp_ConservDF["Maximum identity"] = allSp_ConservDF["Maximum identity"].map(lambda x: x.lstrip("").rstrip("%"))
    allSp_ConservDF["Minimum identity"] = allSp_ConservDF["Minimum identity"].apply(float)  # converting to float
    allSp_ConservDF["Maximum identity"] = allSp_ConservDF["Maximum identity"].apply(float)

    # Apply filters based on minimum and maximum identity thresholds
    allSp_ConservDF_idMin = allSp_ConservDF[allSp_ConservDF["Minimum identity"] >= min_ID]
    if allSp_ConservDF_idMin.empty:  # Check if the DataFrame is empty after filtering
        return allSp_ConservDF_idMin

    allSp_ConservDF_idMax = allSp_ConservDF_idMin[allSp_ConservDF_idMin["Maximum identity"] <= max_ID]
    if allSp_ConservDF_idMax.empty:  # Check if the DataFrame is empty after filtering
        return allSp_ConservDF_idMax

    # Split and process the "Percent of protein sequence matches" column
    seqMatchSplit = allSp_ConservDF_idMax["Percent of protein sequence matches at identity <= 100%"].str.split(" ", n=1, expand=True)
    allSp_ConservDF_idMax["Percent"] = seqMatchSplit[0].apply(float)
    allSp_ConservDF_idMax["Reason"] = seqMatchSplit[1]

    # Filter based on the minimum sequence match percentage
    allSp_ConservDF_final = allSp_ConservDF_idMax[allSp_ConservDF_idMax["Percent"] >= seq_match]

    # Rename columns for clarity
    new_colNames = {
        "Percent of protein sequence matches at identity <= 100%": f"Protein(s) sequence match(es) at Sequence identity threshold {symbol}{ID_threshold}",
        "Minimum identity": "Minimum identity(%)",
        "Maximum identity": "Maximum identity(%)",
    }
    allSp_ConservDF_final = allSp_ConservDF_final.rename(columns=new_colNames, inplace=False)

    # Select and return the final columns of interest
    final_Conserv_df = allSp_ConservDF_final[[
        "Epitope name", "Epitope sequence", "Epitope length", f"Protein(s) sequence match(es) at Sequence identity threshold {symbol}{ID_threshold}", "Minimum identity(%)", "Maximum identity(%)"
    ]]

    return final_Conserv_df
