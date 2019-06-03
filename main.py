from ovrf_functions import get_reading_frames
from ovrf_functions import sort_orfs
from evol_rates import NUCLEOTIDES
from sequence_info import Sequence
from sequence_info import Nucleotide

import argparse


def get_args(parser):
    parser.add_argument(
        'infile',
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'outfile',
        help='Path to the output file'
    )
    parser.add_argument(
        '-pi',
        help='A vector of stationary nucleotide frequencies'
    )
    parser.add_argument(
        '-bias',
        help='List of 6 mutation rates [AC, AG, AT, CG, CT, GT] assuming time-reversibility.',
        default=[1, 1, 1, 1, 1, 1]
    )

    return parser.parse_args()


def main():
    parser = argparse.ArgumentParser(
        description='Simulates and visualizes the evolution of a sequence through a phylogeny'
    )
    args = get_args(parser)

    with open(args.infile, "r") as in_handle:
        # Skip header if the file is a FASTA file
        if in_handle.readline().startswith(">") or in_handle.readline().startswith("#"):
            s = in_handle.readline().strip().upper()


if __name__ == '__main__':
    main()
