import numpy as np
import pandas as pd
import json
from Bio import SeqIO
import re

#Convert JSON format (output of Bepipred 2.0) to pandas DataFrame
def bp2_JsonAnalysis(file):
    """
    Converts a raw JSON file output from Bepipred 2.0 to a pandas DataFrame.

    Parameters:
        file (str): Path to the JSON file.

    Returns:
        pd.DataFrame: DataFrame containing sequence information and prediction results.
    """

    # Load JSON data
    with open(file, 'r') as rawData:
        dataJS = json.load(rawData)
      
    # Extract sequence IDs
    ID_seqs = list(dataJS['antigens'].keys())

    # Define headers based on available information in the JSON file
    headers = ['antigens'] + list(next(iter(dataJS['antigens'].values())).keys())

    # Initialize DataFrame to store results
    final_df = pd.DataFrame(columns=headers)

    # Populate DataFrame with each sequence's prediction data
    for seq in ID_seqs:
        df1 = pd.DataFrame({'antigens': [seq]})
        for header in headers[1:]:
            df2 = pd.DataFrame({header: dataJS['antigens'][seq].get(header, [])},
                               index=np.arange(len(dataJS['antigens'][seq][header])))
            df1 = pd.concat([df1, df2], axis=1)
        final_df = pd.concat([final_df, df1], ignore_index=True)

    # Remove rows with NaN values in the "AA" column
    final_df = final_df.dropna(subset=['AA']).reset_index(drop=True)
    return final_df

# Update DataFrame with a new row of epitope information
def bp2_dfUpdate(df, line_number, InPos, Epitope, new_df, index):
    """
    Updates a DataFrame with new epitope information.

    Parameters:
        df (pd.DataFrame): Original DataFrame with sequence information.
        line_number (int): Row index of sequence information to update.
        InPos (int): Initial position of the epitope.
        Epitope (str): Epitope sequence.
        new_df (pd.DataFrame): DataFrame to append new information.
        index (int): Index for the new entry in `new_df`.

    Returns:
        pd.DataFrame: Updated DataFrame with new epitope data.
    """
    new_row = [
        df['Specie'][line_number],
        df['Protein'][line_number],
        df['ID_Sequence'][line_number],
        InPos,
        df['Position'][line_number],
        Epitope
    ]
    new_df.loc[index] = new_row
    return new_df

# Extract epitope information from Bepipred-2.0 predictions in a DataFrame
def bp2_AntigenEpitopes(dataframe):
    """
    Extracts epitope data from a DataFrame containing Bepipred-2.0 predictions.

    Parameters:
        dataframe (pd.DataFrame): DataFrame with prediction data.

    Returns:
        pd.DataFrame: DataFrame with assembled antigenic regions.
    """
    slice_df = dataframe[['antigens', 'AA', 'PRED']].copy()
    
    # Standardize and map ID sequence information
    idseq = []
    current_id = None
    for line in slice_df['antigens']:
        current_id = line if pd.notna(line) else current_id
        idseq.append(current_id)
    slice_df['antigens'] = idseq

    # Parse species, protein, and sequence ID information
    slice_df[['Protein', 'Specie', 'ID_Sequence']] = slice_df['antigens'].str.extract(r'(\w+?)_(\w+?)_(\w+)')
   
    # Classify residues as antigenic based on prediction score threshold
    slice_df['Classif'] = np.where(slice_df['PRED'] > 0.5, 'Epitope', '-')
     
    # Assign sequential positions to aa residue within each sequence
    slice_df['Position'] = slice_df.groupby('antigens').cumcount() + 1
    
    # Filter and assemble antigenic regions
    antigens_DF = slice_df[slice_df['Classif'] == 'Epitope'].reset_index(drop=True)
    results_df = pd.DataFrame(columns=['Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence'])
   
    index = 0
    epitope_seq = ''
    firstpos = None

    for i in range(len(antigens_DF)):
        pos = antigens_DF['Position'][i]
        if epitope_seq == '':
            firstpos = pos
        epitope_seq += antigens_DF['AA'][i]
        
        if i == len(antigens_DF) - 1 or antigens_DF['Position'][i + 1] != pos + 1:
            results_df = bp2_dfUpdate(antigens_DF, i, firstpos, epitope_seq, results_df, index)
            index += 1
            epitope_seq = ''
    
    return results_df

# Identify epitopes (uppercase segments) in a sequence
def bp3_FindAntigens(sequence):
    """
    Identifies uppercase segments (antigen sequences) in a given sequence.

    Parameters:
        sequence (str): The sequence to search for uppercase segments.

    Returns:
        list: List of dictionaries with antigen sequence and positions.
    """
    antigens = []
    matches = re.finditer(r'[A-Z]+', sequence)

    for match in matches:
        antigen = {
            'Antigen Sequence': match.group(),
            'start_position': match.start() + 1,
            'end_position': match.end()
        }
        antigens.append(antigen)

    return antigens

# Analyze FASTA file (output of Bepipred 3.0) and extract antigen information
def bp3_FastaAnalysis(fastafile):
    """
    Processes a FASTA file output of Bepipred 3.0, identifies uppercase antigen segments, and records them.

    Parameters:
        fastafile (str): Path to the FASTA file.

    Returns:
        pd.DataFrame: DataFrame with information on each antigen segment.
    """
    data = []

    for record in SeqIO.parse(fastafile, "fasta"):
        seq_str = str(record.seq)
        header = record.id
        sp, prot, idSeqNumber = header.split('_', 2)
        antigenSegments = bp3_FindAntigens(seq_str)
        
        for segment in antigenSegments:
            data.append([sp, prot, idSeqNumber, segment['Antigen Sequence'], segment['start_position'], segment['end_position']])

    return pd.DataFrame(data, columns=["Specie", "Protein", "ID_Sequence", "Peptide Sequence", "Initial Position", "Final Position"])


# Apply minimum and maximum length constraints to epitope data
def finalresultsBepipred(results_df, length_min, length_max, BepipredVersion):
    """
    Applies length constraints to predicted epitopes and sets the method version.

    Parameters:
        results_df (pd.DataFrame): DataFrame with epitope predictions.
        length_min (int): Minimum length for epitopes.
        length_max (int): Maximum length for epitopes.
        BepipredVersion (int): Bepipred version (0 for 2.0, 1 for 3.0).

    Returns:
        pd.DataFrame: Filtered DataFrame with epitopes matching length criteria.
    """
    version = 'Bepipred2.0' if BepipredVersion == 0 else 'Bepipred3.0'
    results_df['Method'] = version
    results_df = results_df[['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']]

    if length_min or length_max:
        results_df = results_df[results_df['Peptide Sequence'].str.len().between(length_min or 1, length_max or float('inf'))]
        results_df.reset_index(drop=True, inplace=True)

    return results_df