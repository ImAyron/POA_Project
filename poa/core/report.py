"""Input validation and the POA1 statistics report.

Logic preserved unchanged from the original ``POA1_v1.0.py`` helper functions.
"""
from __future__ import annotations

import math
from datetime import datetime

from Bio import SeqIO


def checkViability(lenght_seq_dict, dataframe, column_epitopes):
    """
    Ensures no epitope is longer than the shortest protein sequence in the dataset,
    which would break the downstream conservancy analysis.

    Raises:
        Exception: If any epitope is longer than the shortest protein sequence.
    """
    for line in range(len(dataframe[column_epitopes])):
        sequence = str(dataframe.iat[line, -1])
        lenght_epitope = len(sequence)
        for seq in lenght_seq_dict:
            if lenght_seq_dict[seq] < lenght_epitope:
                raise Exception(f"ATTENTION: Some of your epitopes (Ex: {dataframe.iat[line, 4]}_{dataframe.iat[line, 5]} ({dataframe.iat[line, 2]}, {dataframe.iat[line, 1]})) have a length greater than the size of the sequence of at least one of the proteins analyzed. Due to the possibility of causing issues in the subsequent analysis (conservation analysis), please reduce the maximum accepted size for the epitopes or consider excluding the protein with the fewest amino acids from the analysis.")


def checkIntegrity(file):
    """
    Checks a FASTA file for invalid characters (e.g. 'x') in protein sequences.

    Returns:
        tuple: (found_invalid: bool, sequence_id: str)
    """
    verif = False
    for seq_record in SeqIO.parse(file, "fasta"):
        ID_sequence = str(seq_record.id).upper()
        sequence = str(seq_record.seq)
        sequence = sequence.lower()
        if "x" in sequence:
            verif = True
            return (verif, ID_sequence)
    return (verif, "")


def writereport(dataframe, path):
    """
    Generates ``Analysis_report.txt`` summarizing the predicted epitopes.
    """
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open(f"{path}/Analysis_report.txt", "w", encoding="utf-8") as report:
        report.write("Analysis Report: POA - Pipeline de Otimização de Antígenos (Antigen Optimization Pipeline)\n\n")
        report.write(dt_string + "\n\n")
        organisms = list(dataframe["Specie"].unique())  # Number of species
        report.write(f"#Species present in the analysis:\n{len(organisms)}\n")
        proteinsInAnalysis = 0  # Number of proteins
        percentGreater = 0
        percentfewer = 0
        Prot_greaterEpit = []
        greater_Epitopes = -(math.inf)
        Prot_fewerEpit = []
        fewer_Epitopes = (math.inf)
        Total_epitopes = len(dataframe["Peptide Sequence"])
        for organism in organisms:
            df_organism = dataframe.loc[(dataframe["Specie"] == organism)]  # per-species df
            organismProteins = list(df_organism["Protein"].unique())
            proteinsInAnalysis += len(organismProteins)
            for protein in organismProteins:
                df_protein = df_organism.loc[(df_organism["Protein"] == protein)]  # per-protein df
                epitopesinProtein = len(df_protein["Peptide Sequence"])
                percent_epitopesinProtein = (epitopesinProtein * 100) / Total_epitopes
                # proteins with fewer epitopes
                if epitopesinProtein < fewer_Epitopes:
                    fewer_Epitopes = len(df_protein["Peptide Sequence"])
                    percentfewer = percent_epitopesinProtein
                    Prot_fewerEpit = []
                    Prot_fewerEpit.append(f"{protein}/{organism}")
                elif epitopesinProtein == fewer_Epitopes:
                    Prot_fewerEpit.append(f"{protein}/{organism}")
                else:
                    pass
                # proteins with more epitopes
                if epitopesinProtein > greater_Epitopes:
                    greater_Epitopes = len(df_protein["Peptide Sequence"])
                    percentGreater = percent_epitopesinProtein
                    Prot_greaterEpit = []
                    Prot_greaterEpit.append(f"{protein}/{organism}")
                elif epitopesinProtein == greater_Epitopes:
                    Prot_greaterEpit.append(f"{protein}/{organism}")
                else:
                    pass
        report.write(f"#Proteins present in the analysis:\n{proteinsInAnalysis}\n")
        report.write("#Epitopes\n\n")
        methods = list(dataframe["Method"].unique())
        for method in methods:
            df_method = dataframe.loc[(dataframe["Method"] == method)]
            epitopesInMethod = len(df_method["Peptide Sequence"])
            report.write(f"{method}:\t{epitopesInMethod}\n")
        report.write(f"\nTotal:\t{Total_epitopes}\n\n")
        report.write("#Proteins with the greatest amount of epitopes\n")
        report.write(f"{Prot_greaterEpit}: {greater_Epitopes} ({percentGreater:.2f}%)\n\n")
        report.write("#Proteins with the least amount of epitopes\n")
        report.write(f"{Prot_fewerEpit}: {fewer_Epitopes} ({percentfewer:.2f}%)\n\n")
        average_length = 0
        for epitope in dataframe["Peptide Sequence"]:
            average_length += len(epitope)
        average_length /= Total_epitopes
        report.write("#Average length of epitopes:\n")
        report.write(f"{int(average_length)} aa.")
