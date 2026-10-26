"""POA1 command-line interface — identical arguments to the original ``POA1_v1.0.py``."""
from __future__ import annotations

import argparse

from ..core.pipeline import run_poa1
from ..logging_conf import setup_logging


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments (at least one prediction method)")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Print this help message.")
    required.add_argument("-b2", help="Bepipred-2.0 prediction in JSON format", type=str, default="")
    required.add_argument("-b3", help="Bepipred-3.0 prediction in FASTA format", type=str, default="")
    optional.add_argument("-bmin", help="Min. Length (MERS) of the predicted Bepipred Epitopes  (default = 0)", type=int, default=0)
    optional.add_argument("-bmax", help="Max. Length (MERS) of the predicted Bepipred Epitopes  (default = 0)", type=int, default=0)
    required.add_argument("-p", help="PAP-IMED prediction in TXT format", type=str, default="")
    optional.add_argument("-pmin", help="Min. Length (MERS) of the predicted PAP-IMED Epitopes  (default = 0)", type=int, default=0)
    optional.add_argument("-pmax", help="Max. Length (MERS) of the predicted PAP-IMED Epitopes  (default = 0)", type=int, default=0)
    required.add_argument("-n", help="NetCTL prediction in HTML format", type=str, default="")
    required.add_argument("-m", help="Directory for MHC-II Binding Predictions in HTML format", type=str, default="")
    optional.add_argument("-mhla", help="HLA-type allele of the predicted MHC-II Binding Predictions Epitopes  (default = DR)", type=str, default="DR")
    optional.add_argument("-mic", help="IC50 threshold - NN_align 2.3 (default = 50), High binding peptides", type=int, default=50)
    required.add_argument("-x", help="Prediction in others web servers (in FASTA format)", type=str, default="")
    optional.add_argument("-xmin", help="Min. Length (MERS) of the results predicted in others web servers  (default = None)", type=int, default=0)
    optional.add_argument("-xmax", help="Max. Length (MERS) of the results predicted in others web servers  (default = None)", type=int, default=0)
    required.add_argument("-d", help="Directory for analysis results", required="True", type=str, default="")
    required.add_argument("-f", help="Polyproteins/proteins used in predictions analysis in fasta format", required="True", type=str, default="")
    optional.add_argument("-e", help="Analysis results in .xlsx format (y/n)", type=str, default="n")
    return parser


def main(argv=None):
    setup_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    print("PrEpiAn: Predicted Epitopes Analyzer v1.0" + "\n")

    result = run_poa1(args)

    print("Data collected and analyzed...\n")
    print(".\n.\n.\n.\n.\n.\n.\n")
    print("Finished")
    return result


if __name__ == "__main__":
    main()
