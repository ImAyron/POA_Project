import pandas as pd
import argparse
import prepianResultsforConservancyAnalysis
import os
import PrEpiAn
import Bio

def checkViability(lenght_seq_dict, dataframe, column_epitopes):
    """
    Verifies the viability of epitopes for conservancy analysis by ensuring that no epitope is longer
    than the shortest protein sequence in the dataset. If an epitope exceeds the length of any protein,
    an exception is raised to prevent errors in downstream analyses.

    Parameters:
        lenght_seq_dict (dict): A dictionary mapping protein IDs to their sequence lengths.
        dataframe (pd.DataFrame): A DataFrame containing the predicted epitopes.
        column_epitopes (str): The name of the column in the DataFrame that contains the epitope sequences.

    Raises:
        Exception: If any epitope is longer than the shortest protein sequence.
    """
    for line in range(len(dataframe[column_epitopes])):
        sequence = str(dataframe.iat[line, -1])
        lenght_epitope = len(sequence)
        for seq in lenght_seq_dict:
            if lenght_seq_dict[seq] < lenght_epitope:
                raise Exception(f"ATTENTION: Some of your epitopes (Ex: {dataframe.iat[line, 4]}_{dataframe.iat[line, 5]} ({dataframe.iat[line, 2]}, {dataframe.iat[line, 1]})) have a length greater than the size of the sequence of at least one of the proteins analyzed. Due to the possibility of causing issues in the subsequent analysis (conservation analysis), please reduce the maximum accepted size for the epitopes or consider excluding the protein with the fewest amino acids from the analysis.")
        
def checkIntegrity(file):
    """
    Checks the integrity of protein sequences in a FASTA file by searching for invalid characters (e.g., 'x').
    If an invalid character is found, the function returns a flag and the ID of the problematic sequence.

    Parameters:
        file (str): Path to the FASTA file containing protein sequences.

    Returns:
        tuple: A tuple containing:
            - verif (bool): True if an invalid character is found, False otherwise.
            - ID_sequence (str): The ID of the sequence containing the invalid character.
    """
    from Bio import SeqIO
    verif = False
    for seq_record in SeqIO.parse(file, "fasta"):
        ID_sequence = str(seq_record.id).upper()
        sequence = str(seq_record.seq)
        sequence = sequence.lower()
        if 'x' in sequence:
            verif = True
            return (verif, ID_sequence)
    return (verif, '')
        
def writereport(dataframe, path):
    """
    Generates a detailed analysis report summarizing the predicted epitopes, including statistics such as:
    - Number of species and proteins analyzed.
    - Total number of epitopes.
    - Proteins with the highest and lowest number of epitopes.
    - Average length of epitopes.

    Parameters:
        dataframe (pd.DataFrame): A DataFrame containing the predicted epitopes and their metadata.
        path (str): Directory where the report will be saved.

    Output:
        A text file named "Analysis_report.txt" containing the analysis summary.
    """
    import math
    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open (f"{path}/Analysis_report.txt", "w") as report:
        report.write("Analysis Report: POA - Pipeline de Otimização de Antígenos (Antigen Optimization Pipeline)\n\n")
        report.write(dt_string+"\n\n")
        organisms = list(dataframe["Specie"].unique())#Number of species
        report.write(f"#Species present in the analysis:\n{len(organisms)}\n")
        proteinsInAnalysis = 0 #Number of proteins
        percentGreater = 0; percentfewer = 0
        Prot_greaterEpit = []; greater_Epitopes = -(math.inf)
        Prot_fewerEpit = []; fewer_Epitopes = (math.inf)
        Total_epitopes = len(dataframe["Peptide Sequence"])
        for organism in organisms:
            df_organism = dataframe.loc[(dataframe["Specie"] == organism)] #dataframe analysis for each species
            organismProteins = list(df_organism["Protein"].unique())
            proteinsInAnalysis += len(organismProteins)
            for protein in organismProteins:
                df_protein = df_organism.loc[(df_organism["Protein"] == protein)] #dataframe analysis by protein 
                epitopesinProtein = len(df_protein["Peptide Sequence"])
                percent_epitopesinProtein = (epitopesinProtein*100)/Total_epitopes
                #proteins with fewer epitopes 
                if epitopesinProtein < fewer_Epitopes:
                    fewer_Epitopes = len(df_protein["Peptide Sequence"])
                    percentfewer = percent_epitopesinProtein
                    Prot_fewerEpit = []
                    Prot_fewerEpit.append(f"{protein}/{organism}")
                elif epitopesinProtein == fewer_Epitopes:
                    Prot_fewerEpit.append(f"{protein}/{organism}")
                else:
                    pass
                #proteins with more epitopes 
                if epitopesinProtein > greater_Epitopes:
                    greater_Epitopes = len(df_protein["Peptide Sequence"])
                    percentGreater = percent_epitopesinProtein
                    Prot_greaterEpit = []
                    Prot_greaterEpit.append(f"{protein}/{organism}")
                elif epitopesinProtein == greater_Epitopes:
                    Prot_greaterEpit.append(f"{protein}/{organism}")
                else:
                    pass
        report.write(f"#Proteins present in the analysis:\n{proteinsInAnalysis}\n")
        report.write("#Epitopes\n\n")
        methods = list(dataframe["Method"].unique())
        for method in methods:
            df_method = dataframe.loc[(dataframe["Method"] == method)]
            epitopesInMethod = len(df_method["Peptide Sequence"])
            report.write(f"{method}:\t{epitopesInMethod}\n")
        report.write(f"\nTotal:\t{Total_epitopes}\n\n")
        report.write("#Proteins with the greatest amount of epitopes\n")
        report.write(f"{Prot_greaterEpit}: {greater_Epitopes} ({percentGreater:.2f}%)\n\n")
        report.write("#Proteins with the least amount of epitopes\n")
        report.write(f"{Prot_fewerEpit}: {fewer_Epitopes} ({percentfewer:.2f}%)\n\n")  
        average_length = 0
        for epitope in dataframe["Peptide Sequence"]:
            average_length += len(epitope)
        average_length /= Total_epitopes
        report.write("#Average length of epitopes:\n")
        report.write(f"{int(average_length)} aa.")
    
def main():
    """
    Main function of the POA1, first stage of POA pipeline. Parses command-line arguments, validates input data, and orchestrates
    the analysis of predicted epitopes. It also generates an analysis report and prepares data for conservancy analysis.

    Command-line arguments:
        -b2: Bepipred-2.0 prediction results in JSON format.
        -b3: Bepipred-3.0 prediction results in FASTA format.
        -bmin: Minimum length of Bepipred epitopes (default = 0).
        -bmax: Maximum length of Bepipred epitopes (default = 0).
        -p: PAP-IMED prediction results in TXT format.
        -pmin: Minimum length of PAP-IMED epitopes (default = 0).
        -pmax: Maximum length of PAP-IMED epitopes (default = 0).
        -n: NetCTL prediction results in HTML format.
        -m: Directory containing MHC-II Binding Predictions in HTML format.
        -mhla: HLA-type allele for MHC-II predictions (default = DR).
        -mic: IC50 threshold for NN_align 2.3 (default = 50).
        -x: Prediction results from other tools in FASTA format.
        -xmin: Minimum length of epitopes from other tools (default = 0).
        -xmax: Maximum length of epitopes from other tools (default = 0).
        -d: Directory to store analysis results.
        -f: FASTA file containing polyproteins/proteins used in predictions.
        -e: Option to export results in .xlsx format (default = 'n').

    Raises:
        Exception: If no prediction method is provided or if invalid characters are found in the protein sequences.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments (at least one prediction method)')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "Print this help message.")
    required.add_argument("-b2", help= "Bepipred-2.0 prediction in JSON format", type=str, default='')
    required.add_argument("-b3", help= "Bepipred-3.0 prediction in FASTA format", type=str, default='')
    optional.add_argument("-bmin", help="Min. Length (MERS) of the predicted Bepipred Epitopes  (default = 0)", type=int, default=0)
    optional.add_argument("-bmax", help="Max. Length (MERS) of the predicted Bepipred Epitopes  (default = 0)", type=int, default=0)
    required.add_argument("-p", help= "PAP-IMED prediction in TXT format", type=str, default='')
    optional.add_argument("-pmin", help="Min. Length (MERS) of the predicted PAP-IMED Epitopes  (default = 0)", type=int, default=0)
    optional.add_argument("-pmax", help="Max. Length (MERS) of the predicted PAP-IMED Epitopes  (default = 0)", type=int, default=0)
    required.add_argument("-n", help= "NetCTL prediction in HTML format", type=str, default='')
    required.add_argument("-m", help= "Directory for MHC-II Binding Predictions in HTML format", type=str, default='')
    optional.add_argument("-mhla", help="HLA-type allele of the predicted MHC-II Binding Predictions Epitopes  (default = DR)", type=str, default = 'DR')
    optional.add_argument("-mic", help="IC50 threshold - NN_align 2.3 (default = 50), High binding peptides", type=int, default = 50) 
    required.add_argument("-x", help= "Prediction in others web servers (in FASTA format)", type=str, default='')
    optional.add_argument("-xmin", help="Min. Length (MERS) of the results predicted in others web servers  (default = None)", type=int, default=0)
    optional.add_argument("-xmax", help="Max. Length (MERS) of the results predicted in others web servers  (default = None)", type=int, default=0)
    required.add_argument("-d", help= "Directory for analysis results",required = "True", type=str, default='')
    required.add_argument("-f", help= "Polyproteins/proteins used in predictions analysis in fasta format", required = "True", type=str, default='')
    optional.add_argument("-e", help= "Analysis results in .xlsx format (y/n)", type=str, default='n')
    
    args = parser.parse_args()
    print("PrEpiAn: Predicted Epitopes Analyzer v1.0" + "\n")
    
    #At least one argument must be selected.
    if (args.b2 == '') and (args.b3 == '') and (args.p == '') and (args.n == '') and (args.m == '') and (args.x == ''):
        raise Exception(f'ERROR: Some the following arguments are required: -n, -m, -b (-b2 or -b3), -p')
    
    #Checking for the presence of 'x' in the sequences 
    check, seq_X = checkIntegrity(args.f)
    if check == True:
        raise Exception(f"ERROR: The sequence {seq_X} contains an invalid character: 'x'")

    #Obtaining the length of protein sequences to verify epitope viability for conservancy analysis 
    from Bio import SeqIO
    lenght_seq_dict = {}
    for seq_record in SeqIO.parse(args.f, "fasta"):
        ID_sequence = str(seq_record.id).upper()
        sequence = str(seq_record.seq)
        sequence = sequence.lower()
        lenght_seq = len(sequence)
        lenght_seq_dict[ID_sequence] = lenght_seq

    #organize prediction results
    results_df = PrEpiAn.runningPrEpiAn(args)

    #Checking if all epitopes are smaller than proteins in analysis
    checkViability(lenght_seq_dict, results_df, "Peptide Sequence")

    #Analysis Report
    writereport(results_df, args.d)

    #Organize data for Epitope Conservancy Analysis
    prepianResultsforConservancyAnalysis.filesforEptConsAnalysis(results_df, args.d)

    print("Data collected and analyzed...\n")
    print(".\n.\n.\n.\n.\n.\n.\n")
    print("Finished")


if __name__ == '__main__':
    main()
