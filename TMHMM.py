def pyTMHMMpredict(file):
    """
    Predicts transmembrane helices using the pyTMHMM tool for protein sequences in a FASTA file.
    The function processes each sequence, extracts its ID, and returns the TMHMM predictions.

    Parameters:
        file (str): Path to the FASTA file containing protein sequences.

    Returns:
        tuple: A tuple containing two lists:
            - id_seqs_list (list): List of protein IDs.
            - prediction_TMHMM_list (list): List of TMHMM predictions for each protein sequence.

    Note:
        The pyTMHMM tool only accepts protein sequences without FASTA headers.
    """
    import pyTMHMM
    from Bio import SeqIO

    id_seqs_list = []
    prediction_TMHMM_list = []

    for seq_record in SeqIO.parse(file, "fasta"):
        id_seq = (seq_record.id).upper() #Extract and standardize the protein ID
        seq = str(seq_record.seq) #Extract the protein sequence
        #Predict transmembrane helices using pyTMHMM     
        annotation = pyTMHMM.predict(seq, compute_posterior=False)
        annotation = annotation.lower() #Standardize the annotation to lowercase
        id_seqs_list.append(id_seq)
        prediction_TMHMM_list.append(annotation)

    return id_seqs_list, prediction_TMHMM_list

def epitTMHMMcaract (seq_TMHMM, initial_pos, lenght):
    """
    Calculates the percentage of amino acids (aa) in an epitope that are located in outer (O),
    transmembrane (M), and inner (I) regions based on TMHMM predictions.

    Parameters:
        seq_TMHMM (str): TMHMM prediction string for a protein sequence.
        initial_pos (int): Starting position of the epitope in the protein sequence.
        lenght (int): Length of the epitope.

    Returns:
        tuple: A tuple containing the percentages of aa in the epitope in:
            - Outer regions (Outporct).
            - Transmembrane regions (TMporct).
            - Inner regions (Insporct).
    """
    #Extract the TMHMM annotation for the epitope region
    epitope_TMHMM = seq_TMHMM[initial_pos:(initial_pos + lenght)]

    #Count the occurrences of each region type
    out_count = epitope_TMHMM.count('o') #Outer regions
    tm_count = epitope_TMHMM.count('m') #Transmembrane regions
    ins_count = epitope_TMHMM.count('i') #Inner regions

    #Calculate the percentage of each region type
    Outporct = out_count / lenght
    TMporct = tm_count / lenght
    Insporct = ins_count / lenght

    #Round the percentages to 4 decimal places
    return (round(Outporct, 4), round(TMporct, 4), round(Insporct, 4))