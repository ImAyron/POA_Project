"""POA2 command-line interface — identical arguments to the original ``POA2_v1.0.py``."""
from __future__ import annotations

import argparse
import sys

from ..core.pipeline import run_poa2
from ..logging_conf import setup_logging


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, description="POA2 - Conservancy Analysis and TMHMM")
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    mutually_exclusive_group = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive_group.add_argument("-g", help="-g: greater than or equal to (>=) in Sequence identity threshold", type=str)
    mutually_exclusive_group.add_argument("-l", help="-l: less-than sign (<) in Sequence identity threshold", type=str)
    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Print this help message.")
    required.add_argument("-d", help="Directory for results of the Epitope Conservancy Analysis in CSV format", type=str, default="", required="True")
    required.add_argument("-t", help="Sequence identity threshold", type=int, default="", required="True")
    optional.add_argument("-r", help="Directory for analysis results", type=str, default="")
    optional.add_argument("-imin", help="Minimum identity (%%) in Conservancy Analysis", type=int, default="60")
    optional.add_argument("-imax", help="Maximum identity (%%) in Conservancy Analysis", type=int, default="100")
    optional.add_argument("-m", help="Percent of protein sequence matches at identity", type=int, default="60")
    required.add_argument("-f", help=" Protein sequence(s) fasta for prediction of transmembrane helices (TMHMM)", required="True")
    optional.add_argument("-rf", help='''Analysis results in fasta format: [0]all epitopes
                                    [1]epitopes in Outside portion(TMHMM)
                                    [2]epitopes in Transmembrane portion(TMHMM)
                                    [3]epitopes in Inside portion(TMHMM)''', type=int, default=None)
    return parser


def main(argv=None):
    setup_logging()
    parser = build_parser()

    # Print help if no arguments are provided
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(argv)

    print("\n" + "POA2 - Analysis Conservancy" + "\n")

    result = run_poa2(args)

    print("Analysis of Results Completed!\n")
    return result


if __name__ == "__main__":
    main()
