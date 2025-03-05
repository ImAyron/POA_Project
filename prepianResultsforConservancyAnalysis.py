def Epitope_EptConsAnalysis(dataframe, path):
    """
    Creates FASTA files for Epitope Conservancy Analysis, organizing epitopes by species.

    Parameters:
        dataframe (pd.DataFrame): DataFrame containing epitope data with columns:
                                  'Specie', 'Protein', 'Method', 'Initial Position', 
                                  'Final Position', 'Peptide Sequence'.
        path (str): Directory path where the FASTA files will be saved.
    """
    #List unique species in the DataFrame
    Specieslist = list(dataframe["Specie"].unique())

    #Process each species
    for sp in Specieslist:
        #Filter epitopes for the current species 
        df_Specie = dataframe.loc[(dataframe['Specie'] == sp)].reset_index(drop = True)

        #Create a FASTA file for the species
        with open (f'{path}/Conservancy Analysis/{sp}_epitopes.fasta', 'w') as epitopes:
            for linha in range(len(df_Specie["Peptide Sequence"])):
                #Create the FASTA header
                id_epitope = f">{df_Specie['Specie'][linha]}_{df_Specie['Protein'][linha]}_{df_Specie['Method'][linha]}_{df_Specie['Initial Position'][linha]}_{df_Specie['Final Position'][linha]}"
                #Write the header and sequence to the file
                epitopes.write(f"{id_epitope}\n")
                epitopes.write(f"{df_Specie['Peptide Sequence'][linha]}\n")
                    
def filesforEptConsAnalysis(EpitopeDataframe, directory):
    """
    Prepares the directory and creates FASTA files for Epitope Conservancy Analysis.

    Parameters:
        EpitopeDataframe (pd.DataFrame): DataFrame containing epitope data.
        directory (str): Directory path where the FASTA files will be saved.
    """
    import os

    #Check if the "Conservancy Analysis" directory exists, and create it if not
    Analysis_directory = False
    for folders, subfolders, files in os.walk(directory):
        if folders == directory:
            if "Conservancy Analysis" in subfolders:
                Analysis_directory = True
        else:
            continue

    if not Analysis_directory:
        dir = f'{directory}/Conservancy Analysis'       
        os.mkdir(dir)

    #Create FASTA files for epitope conservancy analysis
    Epitope_EptConsAnalysis(EpitopeDataframe, directory)