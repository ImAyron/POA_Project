def fasta_epitopes(fasta, lenght_min, lenght_max):
    """
    Extracts and organizes epitopes from a FASTA file, applying optional length constraints.

    Parameters:
        fasta (str): Path to the FASTA file containing epitope sequences.
        lenght_min (int): Minimum length for epitopes. If 0, no minimum length is applied.
        lenght_max (int): Maximum length for epitopes. If 0, no maximum length is applied.

    Returns:
        pd.DataFrame: A DataFrame containing the following columns:
                      - 'Method': Prediction method used.
                      - 'Specie': Species of the protein.
                      - 'Protein': Protein identifier.
                      - 'ID_Sequence': Sequence identifier.
                      - 'Initial Position': Starting position of the epitope.
                      - 'Final Position': Ending position of the epitope.
                      - 'Peptide Sequence': Epitope sequence.
    """
    import pandas as pd

    #Create a new DataFrame to store results
    results_df = pd.DataFrame(columns = ['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence'])

    #Open and parse the FASTA file
    with open (fasta, 'r') as Input:
        index = 0
        for line in Input:
            if not line.isspace():
                line = line.upper().strip()
                if line[0] == '>': #Extract metadata from the header
                    ID_data = line.replace('>', '')
                    ID_data = ID_data.split('_')
                    sp = ID_data[1] #Species
                    prot = ID_data[0] #Protein
                    method = ID_data[2] #Prediction method
                    idSeq = f"{ID_data[3]}_{ID_data[4]}" if ID_data[3] == 'NP' else ID_data[3]  #Sequence ID
                    InitialPos = ID_data[-2] #Initial position
                    FinalPos = ID_data[-1] #Final position
                else: #Extract epitope sequence
                    Epitope = line
                    #Add epitope data to the DataFrame                  
                    results_df.loc[index] = [
                        method, sp, prot, idSeq, InitialPos, FinalPos, Epitope
                    ]
                    index += 1
    #Apply length constraints to epitopes
    if lenght_min == 0 and lenght_max == 0:
        pass #No length filtering
    elif lenght_min == 0:
        results_df = results_df[results_df['Peptide Sequence'].str.len() <= lenght_max]
    elif lenght_max == 0:
        results_df = results_df[results_df['Peptide Sequence'].str.len() >= lenght_min]
    else:
        results_df = results_df[
            (results_df['Peptide Sequence'].str.len() >= lenght_min) & 
            (results_df['Peptide Sequence'].str.len() <= lenght_max)
        ]
        
    results_df = results_df.reset_index(drop = True)
    
    return results_df