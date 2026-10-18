"""Input validation and the POA1 statistics report.

Logic preserved unchanged from the original ``POA1_v1.0.py`` helper functions.
"""
from __future__ import annotations

import math
from datetime import datetime

from Bio import SeqIO


def checkViability(lenght_seq_dict, dataframe, column_epitopes):
    """
    Ensures every epitope fits inside the reference it will later be compared against,
    which would otherwise break the downstream conservancy analysis.

    An epitope is measured against the reference carrying its own ``Protein_Specie`` — comparing
    it against unrelated proteins only produced false refusals. When the reference set holds
    nothing under that ``Protein_Specie``, there is no scoped limit to apply, so the original
    all-references rule stands in: conservancy drops every protein shorter than the epitope from
    that epitope's denominator, and an oversized epitope with no reference of its own would come
    out as ``0.00% (0/0)`` instead of stopping the run.

    Raises:
        ValueError: If any epitope is longer than the reference it will be compared against.
    """
    if not lenght_seq_dict:
        return
    shortest_reference = min(lenght_seq_dict.values())
    for _, row in dataframe.iterrows():
        protein, specie = str(row["Protein"]).upper(), str(row["Specie"]).upper()
        candidates = {key: size for key, size in lenght_seq_dict.items()
                      if key.upper().split("_")[:2] == [protein, specie]}
        reference_id = str(row.get("ID_Sequence", "")).upper()
        exact = f"{protein}_{specie}_{reference_id}"
        if exact in candidates:
            candidates = {exact: candidates[exact]}
        lenght_epitope = len(str(row[column_epitopes]))
        if candidates:
            if lenght_epitope > max(candidates.values()):
                raise ValueError(f"Epitope of {protein}_{specie} has length greater than its matching reference.")
        elif lenght_epitope > shortest_reference:
            raise ValueError(
                f"Epitope of {protein}_{specie} has no matching reference in the -f set and is "
                f"longer ({lenght_epitope} aa) than the shortest reference there "
                f"({shortest_reference} aa). Add a {protein}_{specie}_<ID> reference, or check "
                f"that the predictions and the reference FASTA name this protein the same way."
            )


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
        if Total_epitopes == 0:
            greater_Epitopes = fewer_Epitopes = 0
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
        average_length = average_length / Total_epitopes if Total_epitopes else 0
        report.write("#Average length of epitopes:\n")
        report.write(f"{int(average_length)} aa.")
