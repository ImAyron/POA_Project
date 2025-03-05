#writing the fasta file of the results 
def getfastafile(path,portion_epitopes, dataframe):
        """
    Generates a FASTA file containing epitope sequences based on their membrane topology classification.
    The function writes epitopes to the FASTA file depending on the specified region (outer, transmembrane, or inner).

    Parameters:
        path (str): The directory path where the FASTA file will be saved.
        portion_epitopes (int): Indicates the membrane topology region to include in the FASTA file:
            - 0: All regions.
            - 1: Outer regions.
            - 2: Transmembrane regions.
            - 3: Inner regions.
        dataframe (pd.DataFrame): A DataFrame containing epitope sequences and their membrane topology classification.

    Output:
        A FASTA file named `<path>.fasta` containing the selected epitope sequences.
    """
        with open (f"{path}.fasta", "w") as results_fasta:
            for indice_line in range(len(dataframe["Epitope sequence"])):
                idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
                epitope = dataframe.iloc[indice_line]["Epitope sequence"]
                column = ''
                #Determine which membrane topology region to reference
                if portion_epitopes == 0:#All regions 
                    results_fasta.write(f'>{idt_epitope}')
                    results_fasta.write(f'\n{epitope}\n')
                elif portion_epitopes == 1: #Outer regions
                    column = "Portion_Outside"
                elif portion_epitopes == 2: #Transmembrane regions 
                    column = "Portion_TM"
                else: #Inner regions
                    column = "Portion_Inside"
                if column != '':
                    if dataframe.iloc[indice_line][column] == 1: #Check if the epitope belongs to the specified region
                        results_fasta.write(f'>{idt_epitope}')
                        results_fasta.write(f'\n{epitope}\n')
                 
def main():
    """
    Main function of POA2, second stage of POA pipeline, responsible for performing conservancy analysis and transmembrane helix prediction (TMHMM)
    on epitopes. It processes command-line arguments, organizes the analysis results, and generates output files
    in both Excel and FASTA formats.

    Command-line arguments:
        -g: Greater than or equal to (>=) sequence identity threshold.
        -l: Less than (<) sequence identity threshold.
        -d: Directory containing Epitope Conservancy Analysis results in CSV format.
        -t: Sequence identity threshold (required).
        -r: Directory to store analysis results (optional).
        -imin: Minimum identity percentage in conservancy analysis (default = 60).
        -imax: Maximum identity percentage in conservancy analysis (default = 100).
        -m: Percentage of protein sequence matches at identity (default = 60).
        -f: FASTA file containing protein sequences for TMHMM prediction (required).
        -rf: Option to generate FASTA files for specific membrane topology regions:
            - 0: All epitopes.
            - 1: Epitopes in outer regions.
            - 2: Epitopes in transmembrane regions.
            - 3: Epitopes in inner regions.

    Raises:
        Exception: If required arguments are missing or invalid values are provided.
    """
    import pandas as pd
    import argparse
    import conservancyAnalysis
    import PepiTMHMM
    import sys

    parser = argparse.ArgumentParser(add_help=False, description = 'POA2 - Conservancy Analysis and TMHMM')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    mutually_exclusive_group = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive_group.add_argument("-g", help= "-g: greater than or equal to (>=) in Sequence identity threshold", type=str)
    mutually_exclusive_group.add_argument("-l", help= "-l: less-than sign (<) in Sequence identity threshold", type=str)
    required.add_argument("-h", "--help", action ="help", default=argparse.SUPPRESS, help= "Print this help message.")
    required.add_argument("-d", help= "Directory for results of the Epitope Conservancy Analysis in CSV format", type=str, default='', required = "True")
    required.add_argument("-t", help= "Sequence identity threshold", type=int, default='', required = "True")
    optional.add_argument("-r", help= "Directory for analysis results", type=str, default='')
    optional.add_argument("-imin", help= "Minimum identity (%%) in Conservancy Analysis", type=int, default='60')
    optional.add_argument("-imax", help= "Maximum identity (%%) in Conservancy Analysis", type=int, default='100')
    optional.add_argument("-m", help= "Percent of protein sequence matches at identity", type=int, default='60')
    required.add_argument("-f", help= " Protein sequence(s) fasta for prediction of transmembrane helices (TMHMM)", required = "True")
    optional.add_argument("-rf", help= '''Analysis results in fasta format: [0]all epitopes 
                                    [1]epitopes in Outside portion(TMHMM)
                                    [2]epitopes in Transmembrane portion(TMHMM)
                                    [3]epitopes in Inside portion(TMHMM)''', type=int, default=None)
    # Print help if no arguments are provided
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    print("\n" + 'POA2 - Analysis Conservancy' + "\n")

    # Validate required arguments
    if (args.t == '') or (args.d == ''):
        raise Exception(f'ERROR: The following arguments are required: -t, -d')
    
    # Validate -rf argument
    if args.rf != None:
        if (args.rf < 0) or (args.rf > 3):
            raise Exception(f'ERROR: Argument -rf is invalid')

    # Set output directory path
    if args.r != '':
        path = f"{args.r}/POA2_analysis"
    else:
        path = f"POA2_analysis"

    # Determine the type of conservancy analysis
    if args.g:
        type_symbol = '>='
    else:
        type_symbol = '<'

    #Organize the results of the conservancy analysis 
    ConservancyAnalysis_DF = conservancyAnalysis.EpitConservAnalysis(args.t, type_symbol, args.m, args.imax, args.imin, args.d)   
    
    #Perform TMHMM prediction on epitopes
    POA2_df = PepiTMHMM.tmhmmAnalysis(args, ConservancyAnalysis_DF)

    #Save results to an Excel file
    POA2_df.to_excel(f"{path}_{args.t}.xlsx", sheet_name = f'{type_symbol}{args.t}', startcol=0, index=False)

    #Generate FASTA file if requested
    if (args.rf) != None:
        getfastafile(path, args.rf, POA2_df)

    print("Analysis of Results Completed!\n")

if __name__ == '__main__':
    main()