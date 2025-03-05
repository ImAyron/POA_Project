import pandas as pd

def PAPepitopes(file, length_min=0, length_max=0):
    """
    Extracts Predicted Antigenic Peptides (PAP) from an IMED Server output file,
    filters based on specified peptide length, and organizes results in a DataFrame.

    Parameters:
        file (str): Path to the IMED output file.
        length_min (int, optional): Minimum peptide length for filtering (default is 0, meaning no minimum).
        length_max (int, optional): Maximum peptide length for filtering (default is 0, meaning no maximum).

    Returns:
        pd.DataFrame: DataFrame with columns:
                      'Method', 'Specie', 'Protein', 'ID_Sequence', 
                      'Initial Position', 'Final Position', 'Peptide Sequence'.
    """
    
    # Initialize DataFrame with IMED output structure
    col = ['ID', 'n', 'Start Position', 'Sequence', 'End Position']
    imed_df = pd.DataFrame(columns=col)

    # Parse IMED output file
    with open(file, 'r') as txt:
        data = txt.readlines()
        current_id = ''

        for line in data:
            if line.startswith('>'):
                # Capture identifier information from lines beginning with '>'
                current_id = line.strip().replace('>', '').upper()
            elif line[0] == 'n' or line.isspace():
                # Skip lines that are blank or start with 'n'
                continue
            else:
                # Create a new row entry with ID and peptide sequence information
                row = [current_id] + line.strip().split('\t')
                imed_df.loc[len(imed_df)] = row

    # Organize results into a standardized DataFrame
    results_df = pd.DataFrame(columns=[
        'Method', 'Specie', 'Protein', 'ID_Sequence', 
        'Initial Position', 'Final Position', 'Peptide Sequence'
    ])
 
    # Populate results_df with processed data
    for idx, row in imed_df.iterrows():
        imed_id = row['ID'].split('_')
        specie = imed_id[1] if len(imed_id) > 1 else '-'
        protein = imed_id[0]
        id_seq = "_".join(imed_id[2:4]) if len(imed_id) > 3 and imed_id[2] == "NP" else imed_id[2] if len(imed_id) > 2 else '-'

        # Append the row to results_df with standardized format
        results_df.loc[idx] = [
            'PAP/IMED',
            specie,
            protein,
            id_seq,
            row['Start Position'],
            row['End Position'],
            row['Sequence']
        ]

    # Filter DataFrame by peptide length if specified
    if length_min or length_max:
        results_df = results_df[
            (results_df['Peptide Sequence'].str.len() >= length_min if length_min else True) &
            (results_df['Peptide Sequence'].str.len() <= length_max if length_max else True)
        ].reset_index(drop=True)

    return results_df