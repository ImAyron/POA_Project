"""Pipeline orchestration for POA1 and POA2 — pure logic, reusable by CLI and GUI.

``consolidate_predictions`` is the former ``PrEpiAn.runningPrEpiAn``; ``run_poa1`` / ``run_poa2``
mirror the two original ``main`` functions minus argument parsing. Scientific logic unchanged;
console ``print`` calls became ``logger`` calls.

All functions accept an ``args``-like object (anything with the expected attributes — an
``argparse.Namespace`` from the CLI or a ``SimpleNamespace`` from the GUI).
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

import pandas as pd
from Bio import SeqIO

from ..logging_conf import get_logger
from . import fasta_out, report, topology
from .models import STANDARD_COLUMNS
from .parsers import bepipred, conservancy, mhcii, netctl, others, papimed

logger = get_logger("pipeline")


# --------------------------------------------------------------------------- results containers
@dataclass
class Poa1Result:
    predictions: pd.DataFrame
    report_path: str
    conservancy_dir: str


@dataclass
class Poa2Result:
    results: pd.DataFrame
    xlsx_path: str
    fasta_path: Optional[str]
    type_symbol: str


# --------------------------------------------------------------------------- helpers
def _check_epitope_integrity(dataframe, column_epitopes, column_ID):
    """Raise if any epitope sequence contains an invalid character ('x')."""
    indice = 0
    for line in dataframe[column_epitopes]:
        line = str(line)
        line = line.lower()
        if "x" in line:
            raise Exception(f"ERROR: There is an epitope in sequence {dataframe.iloc[indice][column_ID]} containing an invalid character: 'x'")
        else:
            indice += 1


# --------------------------------------------------------------------------- POA1
def consolidate_predictions(args) -> pd.DataFrame:
    """
    Parse every supplied prediction file and consolidate into one standardized DataFrame.
    (Former ``PrEpiAn.runningPrEpiAn``.)
    """
    import os

    predictions = []

    # Analysis of epitopes predicted by Bepipred
    if (args.b2) != "" or (args.b3) != "":
        if (args.b2) != "" and (args.b3) != "":
            raise Exception("ERROR: Arguments -b2 and -b3 cannot be used simultaneously. Please select results from only one version of Bepipred for analysis at a time.")
        else:
            if (args.b2) != "":
                logger.info("Collecting Bepipred-2.0 epitope data ...")
                bp_data = bepipred.bp2_JsonAnalysis(args.b2)
                Bp_df = bepipred.bp2_AntigenEpitopes(bp_data)
                Bp_df = bepipred.finalresultsBepipred(Bp_df, args.bmin, args.bmax, 0)
            else:
                logger.info("Collecting Bepipred-3.0 epitope data ...")
                Bp_df = bepipred.bp3_FastaAnalysis(args.b3)
                Bp_df = bepipred.finalresultsBepipred(Bp_df, args.bmin, args.bmax, 1)

        _check_epitope_integrity(Bp_df, "Peptide Sequence", "ID_Sequence")
        predictions.append(Bp_df)

    # Analysis of epitopes predicted by PAP (Predicting Antigenic Peptides) - IMED
    if (args.p) != "":
        logger.info("Collecting PAP epitope data ...")
        if (args.pmin) == 0 and (args.pmax) == 0:
            PAP_df = papimed.PAPepitopes(args.p, 0, 0)
        else:
            PAP_df = papimed.PAPepitopes(args.p, args.pmin, args.pmax)
        _check_epitope_integrity(PAP_df, "Peptide Sequence", "ID_Sequence")
        predictions.append(PAP_df)

    # Analysis of epitopes predicted by NetCTL 1.2
    if (args.n) != "":
        logger.info("Collecting NetCTL epitope data ...")
        NetCTL_df = netctl.netctlAntigenEpitopes(args.n)
        _check_epitope_integrity(NetCTL_df, "Peptide Sequence", "ID_Sequence")
        predictions.append(NetCTL_df)

    # Analysis of epitopes predicted by MHCII Binding Prediction (IEDB)
    if (args.m) != "":
        logger.info("Collecting MHCIIBP epitope data ...")
        directory = rf"{args.m}"
        specie, protein, MHCII_files = mhcii.files_map_MHCIIBD(directory)

        columns = ["Method", "Specie", "Protein", "Allele", "ID_Sequence", "Initial Position", "Final Position", "Peptide Sequence"]
        MHCII_df = pd.DataFrame(columns=columns)

        MHCII_df_list = []
        for file in MHCII_files:
            df = mhcii.MHCIIHTMLconverter(file)
            MHCII_df_list.append(df)
        for indice in range(len(MHCII_df_list)):
            df = mhcii.MHCIIAntigenEpitopes(MHCII_df_list[indice], specie[indice], protein[indice], args.mhla, args.mic)
            MHCII_df = pd.concat([MHCII_df, df])
        _check_epitope_integrity(MHCII_df, "Peptide Sequence", "ID_Sequence")
        predictions.append(MHCII_df)

    # Analysis of epitopes predicted by other methods
    if (args.x) != "":
        filename = str(os.path.basename(rf"{args.x}"))
        filename = filename.split(".")
        logger.info("Analyzing data from other predictors ...")
        if (args.xmin) == 0 and (args.xmax) == 0:
            epitopes_df = others.fasta_epitopes(args.x, 0, 0)
        else:
            epitopes_df = others.fasta_epitopes(args.x, args.xmin, args.xmax)
        _check_epitope_integrity(epitopes_df, "Peptide Sequence", "ID_Sequence")
        predictions.append(epitopes_df)

    # Consolidate all predictions into a single DataFrame
    df_predictions = pd.DataFrame(columns=STANDARD_COLUMNS)
    for prediction in predictions:
        df_predictions = pd.concat([df_predictions, prediction], join="inner")

    # Save results to Excel files if requested
    if (args.e).lower() == "y":
        if args.b2 != "" or args.b3 != "":
            Bp_df.to_excel(rf"{args.d}/Bepipred_Epitopes.xlsx", sheet_name="Prediction Results", startcol=0, index=False)
        if args.p != "":
            PAP_df.to_excel(rf"{args.d}/PAP_IMED_Epitopes.xlsx", sheet_name="Prediction Results", startcol=0, index=False)
        if args.n != "":
            NetCTL_df.to_excel(rf"{args.d}/NetCTL_Epitopes.xlsx", sheet_name="Prediction Results", startcol=0, index=False)
        if args.m != "":
            MHCII_df.to_excel(rf"{args.d}/MHCII-BP_Epitopes.xlsx", sheet_name="Prediction Results", startcol=0, index=False)
        if args.x != "":
            epitopes_df.to_excel(rf"{args.d}/{filename[0]}_Epitopes.xlsx", sheet_name="Prediction Results", startcol=0, index=False)

    return df_predictions


def run_poa1(args) -> Poa1Result:
    """
    Full POA1 stage: validate input, consolidate predictions, write the report and the
    per-species FASTA files for the Conservancy Analysis. (Former ``POA1_v1.0.main`` body.)
    """
    # At least one prediction argument must be selected.
    if (args.b2 == "") and (args.b3 == "") and (args.p == "") and (args.n == "") and (args.m == "") and (args.x == ""):
        raise Exception("ERROR: Some the following arguments are required: -n, -m, -b (-b2 or -b3), -p")

    # Ensure the output directory exists (robustness: the pipeline no longer requires the user to
    # pre-create -d; this also covers the optional per-method .xlsx export inside consolidation).
    os.makedirs(args.d, exist_ok=True)

    # Checking for the presence of 'x' in the sequences
    check, seq_X = report.checkIntegrity(args.f)
    if check == True:
        raise Exception(f"ERROR: The sequence {seq_X} contains an invalid character: 'x'")

    # Obtaining the length of protein sequences to verify epitope viability
    lenght_seq_dict = {}
    for seq_record in SeqIO.parse(args.f, "fasta"):
        ID_sequence = str(seq_record.id).upper()
        sequence = str(seq_record.seq).lower()
        lenght_seq_dict[ID_sequence] = len(sequence)

    # Organize prediction results
    results_df = consolidate_predictions(args)

    # Checking if all epitopes are smaller than proteins in analysis
    report.checkViability(lenght_seq_dict, results_df, "Peptide Sequence")

    # Analysis Report
    report.writereport(results_df, args.d)

    # Organize data for Epitope Conservancy Analysis
    fasta_out.filesforEptConsAnalysis(results_df, args.d)

    logger.info("Data collected and analyzed. POA1 finished.")
    return Poa1Result(
        predictions=results_df,
        report_path=f"{args.d}/Analysis_report.txt",
        conservancy_dir=f"{args.d}/Conservancy Analysis",
    )


# --------------------------------------------------------------------------- POA2
def run_poa2(args) -> Poa2Result:
    """
    Full POA2 stage: conservancy filtering + TMHMM topology + xlsx/fasta output.
    (Former ``POA2_v1.0.main`` body.)
    """
    # Validate required arguments
    if (args.t == "") or (args.d == ""):
        raise Exception("ERROR: The following arguments are required: -t, -d")

    # Validate -rf argument
    if args.rf is not None:
        if (args.rf < 0) or (args.rf > 3):
            raise Exception("ERROR: Argument -rf is invalid")

    # Set output directory path
    if args.r != "":
        path = f"{args.r}/POA2_analysis"
    else:
        path = "POA2_analysis"

    # Determine the type of conservancy analysis
    if args.g:
        type_symbol = ">="
    else:
        type_symbol = "<"

    # Organize the results of the conservancy analysis
    ConservancyAnalysis_DF = conservancy.EpitConservAnalysis(args.t, type_symbol, args.m, args.imax, args.imin, args.d)

    # Perform TMHMM prediction on epitopes
    POA2_df = topology.tmhmmAnalysis(args, ConservancyAnalysis_DF)

    # Save results to an Excel file (create the output directory if needed)
    xlsx_path = f"{path}_{args.t}.xlsx"
    os.makedirs(os.path.dirname(xlsx_path) or ".", exist_ok=True)
    POA2_df.to_excel(xlsx_path, sheet_name=f"{type_symbol}{args.t}", startcol=0, index=False)

    # Generate FASTA file if requested
    fasta_path = None
    if (args.rf) is not None:
        fasta_out.getfastafile(path, args.rf, POA2_df)
        fasta_path = f"{path}.fasta"

    logger.info("Analysis of Results Completed!")
    return Poa2Result(results=POA2_df, xlsx_path=xlsx_path, fasta_path=fasta_path, type_symbol=type_symbol)
