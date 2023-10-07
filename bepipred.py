#Convert json format (Bepipred's output) to pandas DataFrame 
def BpJsonAnalysis(file):
    import numpy as np
    import pandas as pd
    import json

    #Get json file information
    with open(file, 'r') as rawData:
        dataJS = json.load(rawData)
      
    #Obtain the species names
    org = []
    for x in dataJS['antigens']:
        org.append(x)
      
    #Get header information from the json file
    data = ['antigens']
    for specie in org:
        for info in dataJS['antigens'][specie]:
            if info not in data:
                data.append(info)
    
    #Create a new DataFrame for results
    final_df = pd.DataFrame()
    final_df = pd.DataFrame(columns=data)

    #Add info from each prediction (the file can contain results from multiple predictions) in the new DataFrame
    for specie in org:
        #Create and set a Dataframe for each prediction data (prediction for a single sequence)
        for info in data:
            if info == 'antigens':
                df = pd.DataFrame({f'{info}':[specie]}) 
            else:
                df1 = pd.DataFrame({f'{info}':dataJS['antigens'][specie][info]}, index = np.arange(len(dataJS['antigens'][specie][info]))) 
                df = pd.concat([df,df1], ignore_index=False, axis=1)
        #Merge the two Dataframes (Dataframe for prediction data of single sequence and global Dataframe)
        final_df = pd.merge(final_df, df, how = 'outer')
    
    #Remove any rows with NaN values in antigen columns
    mask = final_df['AA'].isna()
    final_df = final_df.loc[~mask]
    final_df.reset_index(drop=True, inplace=True)

    return(final_df)

#Dataframe update with new lines
def dfUpdate(df, line_number, InPos, Epitope, new_df, index):
    new_line = [    
            df['Specie'][line_number],
            df['Protein'][line_number],
            df['ID_Sequence'][line_number],
            InPos,
            df['Position'][line_number],
            Epitope
        ]
    new_df.loc[index] = new_line
    return new_df
#Get predicted epitope information from BepiPred-2.0 Server from DataFrame pandas
def bpAntigenEpitopes(dataframe, lenght_min, lenght_max):
    import pandas as pd
    #Extract data of interest on DataFrame
    raw_df = dataframe
    slice_df = raw_df[['antigens', 'AA', 'PRED']]
    #Organize sequence information
    idseq = []
    for line in raw_df['antigens']:
        line = str(line) #data standardization
        line.lower() #data standardization
        if line != 'nan':#Get data from sequences
            etiqueta = line
            idseq.append(line) 
        else: #Map positions with Nan and replace them with sequence data
            idseq.append(etiqueta)
    slice_df = slice_df.assign(antigens = idseq) #Update the column with collected data
    
    #Get information on species, proteins and ID sequence via the "antigens" column.
    prot = [];  sp = []; idSeqNumber = []
    for line in slice_df['antigens']:
        dataID = line.split('_')
        prot.append(dataID[0]) 
        sp.append(dataID[1])
        if len(dataID) >= 4:
            idSeqNumber.append(dataID[2] + "_" + dataID[3])
        else:
            idSeqNumber.append(dataID[2])
    
    #Add new columns in the Dataframe
    slice_df['Specie'] = sp
    slice_df['Protein'] = prot
    slice_df['ID_Sequence'] = idSeqNumber
    df = slice_df
    
    #Remove the "antigens" column
    df = df.drop(columns='antigens')
    #Selecting the score for antigenic amino acid prediction
    antigenic_list = []
    for line in range(len(df['PRED'])):
        if df['PRED'][line] > 0.5:
            antigenic_list.append('Epitope')
        else:
            antigenic_list.append('-')
    
    #Add a new column to characterize amino acid residues that are antigenic
    df['Classif'] = antigenic_list

    #Adding a column for amino acid residue position data
    ref = df['Specie'][0]
    pos = []; cont = 0 
    for line in df['Specie']:
        if line == ref:
            cont += 1
            pos.append(cont)
        else:
            ref = line
            pos.append(1)
            cont = 1
    df['Position'] = pos

    #Rearranging the position of columns in DataFrame

    df = df[['Specie', 'Protein', 'ID_Sequence', 'Position', 'AA', 'PRED', 'Classif']]

    #DataFrame of immunogenic residues
    antigens_DF = df['Classif'] == 'Epitope'
    antigens_DF = df[antigens_DF]
    antigens_DF = antigens_DF.reset_index(drop = True)

    #Assembly of antigenic regions (based on amino acid position)
    #New DataFrame for the results
    results_df = pd.DataFrame()
    col = ['Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
    results_df = pd.DataFrame(columns=col)
    #Assembling the epitope and retrieving localization information from the single amino acid residues
    index = 0
    epit = ''
    for line in range(len(antigens_DF['Position'])):
        x = antigens_DF['Position'][line]
        if epit == '':
            firstpos = int(x)
        #Case of the amino acid residues in the last row of the Dataframe
        if line == ((len(antigens_DF['Position']))-1):
            if x - 1 == antigens_DF['Position'][line - 1]: #Looking at the position of the amino acid residue from the previous line
                epit += antigens_DF['AA'][line]
                results_df = dfUpdate(antigens_DF, line, firstpos, epit, results_df, index)
                index += 1; epit = ''
            else:
                epit = antigens_DF['AA'][line]
                results_df = dfUpdate(antigens_DF, line, antigens_DF['Position'][line], epit, results_df, index)
                index += 1; epit = ''

        #Case of the amino acid residues that aren't in the last row
        else:
            if x + 1 == antigens_DF['Position'][line + 1]: #Looking at the position of the amino acid residue in the next line
                epit += antigens_DF['AA'][line]
            elif x - 1 == antigens_DF['Position'][line - 1]: #Looking at the position of the amino acid residue from the previous line
                epit += antigens_DF['AA'][line]
                results_df = dfUpdate(antigens_DF, line, firstpos, epit, results_df, index)
                index += 1; epit = ''
            else:
                epit = antigens_DF['AA'][line]
                results_df = dfUpdate(antigens_DF, line, antigens_DF['Position'][line], epit, results_df, index)
                index += 1; epit = ''
    #Building the final DataFrame
    results_df = results_df.assign(Method='Bepipred')
    results_df = results_df[['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']]
    
    #To set the minimum and maximum length of epitopes
    if lenght_min == 0 and lenght_max == 0:
        results_df = results_df
    elif lenght_min == 0:
        results_df = results_df[results_df['Peptide Sequence'].str.len() <= lenght_max]
        results_df = results_df.reset_index(drop = True)
    elif lenght_max == 0:
        results_df = results_df[results_df['Peptide Sequence'].str.len() >= lenght_min]
        results_df = results_df.reset_index(drop = True)
    else:
        results_df = results_df[results_df['Peptide Sequence'].str.len() >= lenght_min]
        results_df = results_df[results_df['Peptide Sequence'].str.len() <= lenght_max]
        results_df = results_df.reset_index(drop = True)
    return results_df
