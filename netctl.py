import pandas as pd
import html

def netctlAntigenEpitopes(file):
    """
    Extracts predicted epitopes from NetCTL-1.2 Server results in HTML format, 
    creating a DataFrame with standardized columns.
    
    Parameters:
        file (str): Path to the HTML file containing NetCTL prediction results.
        
    Returns:
        pd.DataFrame: DataFrame with columns:
                      'Method', 'Specie', 'Protein', 'ID_Sequence', 
                      'Initial Position', 'Final Position', 'Peptide Sequence'.
    """
    #Define column names according to NetCTL result table and initialize DataFrame
    columns = [
        'Residue_number', 'ID', 'Protein_identifier', 'pep', 'Peptide_sequence', 
        'aff', 'Predicted_MHC_binding_affinity', 'aff_rescale', 'Rescale_binding_affinity',
        'cle', 'C_terminal_cleavage_affinity', 'tap', 'TAP_transport_efficiency',
        'COMB', 'Prediction_score', 'Identified_MHC_ligands'
    ]
    netctl_results_df = pd.DataFrame(columns=columns)

    #Read and parse HTML file contents
    with open(file, 'r') as raw_data:
        netctl_results = [html.unescape(line.strip()) for line in raw_data if not line.isspace()]

    #Fill DataFrame with parsed NetCTL results
    for i, line in enumerate(netctl_results):
        if line and line[0].isdigit():
            line = ' '.join(line.split())  #Remove extra spaces
            new_df_line = line.split(' ')
            
            #If last element is not 'E', add a placeholder
            if new_df_line[-1] != 'E':
                new_df_line.append('-')
            
            netctl_results_df.loc[i] = new_df_line

    #Select predicted epitopes identified by NetCTL
    netctl_epitopes = netctl_results_df[netctl_results_df['Identified_MHC_ligands'] == '<-E'].reset_index(drop=True)

    #Select and rename relevant columns
    netctl_epitopes_slice = netctl_epitopes[['Protein_identifier', 'Residue_number', 'Peptide_sequence']].copy()
    netctl_epitopes_slice.rename(columns={'Residue_number': 'Initial Position', 'Peptide_sequence': 'Peptide Sequence'}, inplace=True)
 
    #Add standard columns
    netctl_epitopes_slice['Method'] = 'NetCTL'
    netctl_epitopes_slice['ID_Sequence'] = '-'

    #Split 'Protein_identifier' to extract 'Specie' and 'Protein'
    netctl_epitopes_slice['Specie'] = netctl_epitopes_slice['Protein_identifier'].apply(lambda x: x.split('_')[1])
    netctl_epitopes_slice['Protein'] = netctl_epitopes_slice['Protein_identifier'].apply(lambda x: x.split('_')[0])

    #Calculate 'Final Position' based on 'Initial Position' and peptide length
    netctl_epitopes_slice['Final Position'] = netctl_epitopes_slice.apply(
        lambda row: int(row['Initial Position']) + len(row['Peptide Sequence']) - 1, axis = 1
    )
    
    #Organize columns in final DataFrame
    results_df = netctl_epitopes_slice[['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']]
    
    return results_df