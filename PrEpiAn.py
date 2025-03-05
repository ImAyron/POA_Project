#Verification of the viability of epitopes to be analyzed by conservancy analysis 
def checkIntegrity(dataframe, column_epitopes, column_ID):
    """
    Verifies the integrity of epitope sequences in a DataFrame by checking for invalid characters (e.g., 'x').
    If an invalid character is found, an exception is raised.

    Parameters:
        dataframe (pd.DataFrame): DataFrame containing the epitope sequences.
        column_epitopes (str): Name of the column containing the epitope sequences.
        column_ID (str): Name of the column containing the sequence IDs.

    Raises:
        Exception: If an invalid character is found in any epitope sequence.
    """
    indice = 0
    for line in dataframe[column_epitopes]:
        line = str(line)
        line = line.lower()
        if 'x' in line:
            raise Exception(f"ERROR: There is an epitope in sequence {dataframe.iloc[indice][column_ID]} containing an invalid character: 'x'")
        else:
            indice += 1

def runningPrEpiAn(args):
    """
    Orchestrates the analysis of epitope predictions from multiple tools (Bepipred, PAP-IMED, NetCTL, MHC-II Binding, and others).
    It processes input files, validates epitope sequences, and consolidates the results into a single DataFrame.

    Parameters:
        args (argparse.Namespace): Command-line arguments containing input file paths and options.

    Returns:
        pd.DataFrame: A consolidated DataFrame containing epitope predictions from all tools.

    Raises:
        Exception: If invalid arguments are provided or if epitope sequences contain invalid characters.
    """
    import pandas as pd
    import mhcii
    import netctl
    import papImed
    import bepipred
    import othersPred
    import os

    predictions = []

    #Analysis of epitopes predicted by Bepipred
    if (args.b2) != '' or (args.b3) != '':
        if (args.b2) != '' and (args.b3) != '':
            raise Exception(f"ERROR: Arguments -b2 and -b3 cannot be used simultaneously. Please select results from only one version of Bepipred for analysis at a time.")
        else:
            if (args.b2) != '':
                print("Collecting Bepipred-2.0 epitope data ...\n")
                bp_data = bepipred.bp2_JsonAnalysis(args.b2)
                Bp_df = bepipred.bp2_AntigenEpitopes(bp_data)
                Bp_df = bepipred.finalresultsBepipred(Bp_df, args.bmin, args.bmax, 0) 
            else:
                print("Collecting Bepipred-3.0 epitope data ...\n")
                Bp_df = bepipred.bp3_FastaAnalysis(args.b3)
                Bp_df = bepipred.finalresultsBepipred(Bp_df, args.bmin, args.bmax, 1) 

        checkIntegrity(Bp_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(Bp_df)
        
    #Analysis of epitopes predicted by PAP (Predicting Antigenic Peptides) - IMED
    if (args.p) != '':
        print("Collecting PAP epitope data ...\n")
        if (args.pmin) == 0 and (args.pmax) == 0:
            PAP_df = papImed.PAPepitopes(args.p, 0, 0)
        else:
            PAP_df = papImed.PAPepitopes(args.p, args.pmin, args.pmax)
        checkIntegrity(PAP_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(PAP_df)
            
    #Analysis of epitopes predicted by NetCTL 1.2
    if (args.n) != '':
        print("Collecting NetCTL epitope data ...\n")
        NetCTL_df = netctl.netctlAntigenEpitopes(args.n)
        checkIntegrity(NetCTL_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(NetCTL_df)
        
    #Analysis of epitopes predicted by MHCII Binding Prediction (IEDB)
    if (args.m) != '':
        print("Collecting MHCIIBP epitope data ...\n")
        directory = rf"{args.m}"
        specie, protein, MHCII_files = mhcii.files_map_MHCIIBD(directory)
        
        MHCII_df = pd.DataFrame()
        columns = ['Method', 'Specie', 'Protein', 'Allele', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
        MHCII_df = pd.DataFrame(columns=columns)

        MHCII_df_list = []
        for file in MHCII_files:
            df = mhcii.MHCIIHTMLconverter(file)
            MHCII_df_list.append(df)
        for indice in range(len(MHCII_df_list)):
            df = mhcii.MHCIIAntigenEpitopes(MHCII_df_list[indice], specie[indice], protein[indice], args.mhla, args.mic)
            MHCII_df = pd.concat([MHCII_df, df])
        checkIntegrity(MHCII_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(MHCII_df)

    #Analysis of epitopes predicted by other methods
    if (args.x) != '':
        filename = str(os.path.basename(rf"{args.x}"))
        filename = filename.split('.')
        print("Analyzing data from other predictors ...\n")
        if (args.xmin) == 0 and (args.xmax) == 0:
            epitopes_df = othersPred.fasta_epitopes(args.x, 0, 0)
        else:
            epitopes_df = othersPred.fasta_epitopes(args.x, args.xmin, args.xmax)
        checkIntegrity(epitopes_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(epitopes_df)
    
    #Consolidate all predictions into a single DataFrame
    df_predictions = pd.DataFrame()
    columns = ['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
    df_predictions = pd.DataFrame(columns=columns)
    for prediction in predictions:
        df_predictions = pd.concat([df_predictions, prediction], join="inner")
    
    #Save results to Excel files if requested
    if (args.e).lower() == 'y':
        if args.b2 != '' or args.b3 != '':
            Bp_df.to_excel(rf'{args.d}/Bepipred_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.p != '':
            PAP_df.to_excel(rf'{args.d}/PAP_IMED_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.n != '':
            NetCTL_df.to_excel(rf'{args.d}/NetCTL_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.m != '':
            MHCII_df.to_excel(rf"{args.d}/MHCII-BP_Epitopes.xlsx", sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.x != '':
            epitopes_df.to_excel(rf'{args.d}/{filename[0]}_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)

    return(df_predictions)