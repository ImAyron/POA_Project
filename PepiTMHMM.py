def notEmptyValidate(file):
    """
    Validates whether a file is empty by checking its size.

    Parameters:
        file (str): Path to the file to be validated.

    Returns:
        bool: True if the file is empty, False otherwise.
    """
    import os
    result = os.stat(file).st_size==0
    return result

#Check if epitope position values are possible 
def rangeVerify(valor, lenghtProt):
    """
    Verifies if a given position value is within the valid range for a protein sequence.

    Parameters:
        valor (int or str): The position value to be validated.
        lenghtProt (int): The length of the protein sequence.

    Returns:
        bool: True if the position is valid, False otherwise.
    """
    if int(valor) == 0 or int(valor) < 0:
        return False
    elif int(valor) > int(lenghtProt):
        return False
    else:
        return True

def tmhmmAnalysis(args, dataframe):
    """
    Performs transmembrane helix prediction (TMHMM) on epitopes and updates the DataFrame with membrane topology classifications.

    Parameters:
        args (argparse.Namespace): Command-line arguments containing input file paths and options.
        dataframe (pd.DataFrame): DataFrame containing epitope sequences and metadata.

    Returns:
        pd.DataFrame: Updated DataFrame with membrane topology classifications for each epitope.

    Raises:
        Exception: If the protein FASTA file is empty or if epitope positions are invalid.
        Warning: If an epitope sequence is not found in the provided protein sequences.
    """
    import TMHMM
    import warnings
    import pandas as pd
    import Bio
    from Bio import SeqIO

    print("TMHMM - Epitope Analysis v1.0" + "\n")

    #Check if the protein FASTA file is empty
    if notEmptyValidate(args.f):
        raise Exception("Failed because the protein FASTA file is empty.")
        
    #Print protein sequences for verification
    for seq_record in SeqIO.parse(args.f, "fasta"):
        print("Completed Sequence Analysis:"+ seq_record.id)       

    #Add columns for membrane topology classification to the DataFrame   
    dataframe = dataframe.assign(Portion_Outside = "-")
    dataframe = dataframe.assign(Portion_TM = "-")
    dataframe = dataframe.assign(Portion_Inside = "-")

    #Perform TMHMM prediction on protein sequences
    ids_TMHMM_list, seqs_TMHMM_list = TMHMM.pyTMHMMpredict(args.f)

    #Lists to store membrane topology classification results
    Out_portion_list = []; TM_portion_list = []; Ins_portion_list = []

    #Analyze each epitope in the DataFrame
    for indice_line in range(len(dataframe["Epitope sequence"])):
        #Extract epitope information 
        idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
        idt = idt_epitope.upper()
        idt = idt.split('_')
        if len(idt) == 5:
            EpitopeVirus = idt[0]
            init_pos = idt[-2]
            fin_pos = idt[-1]
        else:
            raise Exception(f"Failed because epitope ID {idt_epitope} is not correctly formatted. Found {len(idt)} parts instead of 5.")
        
        # Extract and standardize the epitope sequence     
        epitope = dataframe.iloc[indice_line]["Epitope sequence"]
        epitope = epitope.lower()
        
        #Search for matching protein sequences 
        for seq_record in SeqIO.parse(args.f, "fasta"):
            seq_polyprot = str(seq_record.seq).lower()
            seq_polyprot_id = seq_record.id.upper()

            #Check if the protein sequence belongs to the same virus as the epitope
            if EpitopeVirus not in seq_polyprot_id:
                continue

            #Ensure the epitope sequence is present in the protein sequence
            if epitope not in seq_polyprot:
                warnings.warn(f"Epitope {idt_epitope} is not found in the protein(s) FASTA file. Please check your sequences.")
                continue

            # Validate epitope position values
            start_pos_verif = rangeVerify(init_pos, len(seq_record))
            end_pos_verif = rangeVerify(fin_pos, len(seq_record))
            if (start_pos_verif != True) or (end_pos_verif != True):
                raise Exception(f"Failed because epitope positions in {idt_epitope} are invalid.")
                    
            #Match epitope with TMHMM prediction results        
            for x in range(len(seqs_TMHMM_list)):
                if seq_polyprot_id == ids_TMHMM_list[x]:
                    init = seq_polyprot.find(epitope)
                    if init != -1:
                        epit_lenght = len(epitope)
                        #Calculate the percentage of residues in each membrane topology classification 
                        Out_portion, TM_portion, Ins_portion = TMHMM.epitTMHMMcaract (seqs_TMHMM_list[x], init, epit_lenght)
                        Out_portion_list.append(Out_portion)
                        TM_portion_list.append(TM_portion)
                        Ins_portion_list.append(Ins_portion)
    
    #Update DataFrame with membrane topology classifications
    dataframe["Portion_Outside"] = Out_portion_list
    dataframe["Portion_TM"] = TM_portion_list
    dataframe["Portion_Inside"] = Ins_portion_list
    
    return (dataframe)