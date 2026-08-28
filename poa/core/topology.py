"""Membrane-topology characterization of epitopes via pyTMHMM.

Combines the original ``TMHMM.py`` (pyTMHMM wrapper + per-epitope region math) and
``PepiTMHMM.py`` (epitope↔protein matching). Scientific logic preserved unchanged.

``pyTMHMM`` is imported lazily inside :func:`pyTMHMMpredict` so this module can be imported
(for the rest of the pipeline / tests) even when pyTMHMM is not installed.
"""
from __future__ import annotations

import os
import warnings

from Bio import SeqIO

from ..logging_conf import get_logger

logger = get_logger("topology")


# --------------------------------------------------------------------------- pyTMHMM wrapper
def pyTMHMMpredict(file):
    """
    Predicts transmembrane helices for each sequence in a FASTA file using pyTMHMM.

    Returns:
        tuple(list, list): protein IDs (upper-cased) and lowercase o/m/i annotation strings.
    """
    import pyTMHMM  # lazy import: only required for POA2 topology

    id_seqs_list = []
    prediction_TMHMM_list = []

    for seq_record in SeqIO.parse(file, "fasta"):
        id_seq = (seq_record.id).upper()  # Extract and standardize the protein ID
        seq = str(seq_record.seq)  # Extract the protein sequence
        # Predict transmembrane helices using pyTMHMM
        annotation = pyTMHMM.predict(seq, compute_posterior=False)
        annotation = annotation.lower()  # Standardize the annotation to lowercase
        id_seqs_list.append(id_seq)
        prediction_TMHMM_list.append(annotation)

    return id_seqs_list, prediction_TMHMM_list


def epitTMHMMcaract(seq_TMHMM, initial_pos, lenght):
    """
    Fraction of an epitope's residues in outer (o), transmembrane (m) and inner (i) regions.

    Returns:
        tuple(float, float, float): (outer, transmembrane, inner) rounded to 4 decimals.
    """
    # Extract the TMHMM annotation for the epitope region
    epitope_TMHMM = seq_TMHMM[initial_pos:(initial_pos + lenght)]

    # Count the occurrences of each region type
    out_count = epitope_TMHMM.count("o")  # Outer regions
    tm_count = epitope_TMHMM.count("m")  # Transmembrane regions
    ins_count = epitope_TMHMM.count("i")  # Inner regions

    # Calculate the percentage of each region type
    Outporct = out_count / lenght
    TMporct = tm_count / lenght
    Insporct = ins_count / lenght

    # Round the percentages to 4 decimal places
    return (round(Outporct, 4), round(TMporct, 4), round(Insporct, 4))


# --------------------------------------------------------------------------- epitope matching
def notEmptyValidate(file):
    """Return True if ``file`` is empty (size 0)."""
    return os.stat(file).st_size == 0


def rangeVerify(valor, lenghtProt):
    """Verify a position value is within [1, lenghtProt]."""
    if int(valor) == 0 or int(valor) < 0:
        return False
    elif int(valor) > int(lenghtProt):
        return False
    else:
        return True


def sameVirus(epitope_virus, protein_id):
    """
    Whether a ``-f`` record belongs to the same organism as an epitope.

    The epitope name is ``Specie_Protein_Method_Init_Final``, so its first field is the species.
    Protein headers follow ``Protein_Specie_ID``, so the comparison is against their *species
    field*; a plain substring test matches too much — ``DENV1`` is contained in both
    ``DENV1_NS1`` and ``E_DENV1_REF``, which made a single epitope match two records.
    Headers with no ``_`` have no species field and keep the original substring behaviour.
    """
    parts = str(protein_id).upper().split("_")
    if len(parts) > 1:
        return parts[1] == str(epitope_virus).upper()
    return str(epitope_virus).upper() in str(protein_id).upper()


def tmhmmAnalysis(args, dataframe):
    """
    Runs TMHMM on the proteins and annotates each epitope with its membrane-topology fractions.

    Returns:
        pd.DataFrame: input dataframe with Portion_Outside / Portion_TM / Portion_Inside filled.
    """
    logger.info("TMHMM - Epitope Analysis v1.0")

    # Check if the protein FASTA file is empty
    if notEmptyValidate(args.f):
        raise Exception("Failed because the protein FASTA file is empty.")

    # Print protein sequences for verification
    for seq_record in SeqIO.parse(args.f, "fasta"):
        logger.info("Completed Sequence Analysis: %s", seq_record.id)

    # Add columns for membrane topology classification to the DataFrame
    dataframe = dataframe.assign(Portion_Outside="-")
    dataframe = dataframe.assign(Portion_TM="-")
    dataframe = dataframe.assign(Portion_Inside="-")

    # Perform TMHMM prediction on protein sequences
    ids_TMHMM_list, seqs_TMHMM_list = pyTMHMMpredict(args.f)
    annotation_by_id = dict(zip(ids_TMHMM_list, seqs_TMHMM_list))
    records = [(str(rec.id).upper(), str(rec.seq).lower(), len(rec))
               for rec in SeqIO.parse(args.f, "fasta")]

    # Lists to store membrane topology classification results — FIX: exactly ONE entry per
    # DataFrame row. They used to be appended once per *matching protein record*, so an epitope
    # matching two records (or none) desynchronised the lists from the frame: pandas then either
    # rejected the assignment ("Length of values does not match length of index") or, when the
    # extra and missing matches happened to cancel out, silently wrote each row's topology onto a
    # different epitope.
    Out_portion_list = []
    TM_portion_list = []
    Ins_portion_list = []

    # Analyze each epitope in the DataFrame
    for indice_line in range(len(dataframe["Epitope sequence"])):
        # Extract epitope information
        idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
        idt = idt_epitope.upper()
        idt = idt.split("_")
        if len(idt) == 5:
            EpitopeVirus = idt[0]
            init_pos = idt[-2]
            fin_pos = idt[-1]
        else:
            raise Exception(f"Failed because epitope ID {idt_epitope} is not correctly formatted. Found {len(idt)} parts instead of 5.")

        # Extract and standardize the epitope sequence
        epitope = dataframe.iloc[indice_line]["Epitope sequence"]
        epitope = epitope.lower()

        # Search for matching protein sequences
        matches = []
        for seq_polyprot_id, seq_polyprot, rec_lenght in records:
            # Check if the protein sequence belongs to the same virus as the epitope
            if not sameVirus(EpitopeVirus, seq_polyprot_id):
                continue

            # Ensure the epitope sequence is present in the protein sequence
            if epitope not in seq_polyprot:
                warnings.warn(f"Epitope {idt_epitope} is not found in the protein(s) FASTA file. Please check your sequences.")
                continue

            # Validate epitope position values
            start_pos_verif = rangeVerify(init_pos, rec_lenght)
            end_pos_verif = rangeVerify(fin_pos, rec_lenght)
            if (start_pos_verif != True) or (end_pos_verif != True):
                raise Exception(f"Failed because epitope positions in {idt_epitope} are invalid.")

            # Match epitope with TMHMM prediction results
            if seq_polyprot_id in annotation_by_id:
                matches.append((seq_polyprot_id, seq_polyprot))

        if not matches:
            # No protein carries this epitope: keep the '-' placeholder for THIS row only.
            warnings.warn(f"Epitope {idt_epitope} matched no protein in the FASTA file; its membrane topology is reported as '-'.")
            Out_portion_list.append("-")
            TM_portion_list.append("-")
            Ins_portion_list.append("-")
            continue

        if len(matches) > 1:
            warnings.warn(
                f"Epitope {idt_epitope} matched {len(matches)} protein records "
                f"({', '.join(m[0] for m in matches)}); using {matches[0][0]}. Check for duplicated "
                "sequences in the -f FASTA file."
            )

        seq_polyprot_id, seq_polyprot = matches[0]
        init = seq_polyprot.find(epitope)
        epit_lenght = len(epitope)
        # Calculate the percentage of residues in each membrane topology classification
        Out_portion, TM_portion, Ins_portion = epitTMHMMcaract(annotation_by_id[seq_polyprot_id], init, epit_lenght)
        Out_portion_list.append(Out_portion)
        TM_portion_list.append(TM_portion)
        Ins_portion_list.append(Ins_portion)

    # Update DataFrame with membrane topology classifications
    dataframe["Portion_Outside"] = Out_portion_list
    dataframe["Portion_TM"] = TM_portion_list
    dataframe["Portion_Inside"] = Ins_portion_list

    return dataframe
