from ovrf_functions import get_reading_frames
from ovrf_functions import sort_orfs
from evol_rates import NUCLEOTIDES
from sequence_info import Sequence
from sequence_info import Nucleotide

import argparse


def get_args(parser):
    parser.add_argument(
        'sequence',
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'tree',
        help='Path to file containing phylogenetic tree in Newick format.'
    )
    parser.add_argument(
        'rates',
        help='File containing substitution biases to calculate substitution rates'
    )
    parser.add_argument(
        'outfile',
        help='Path to the alignment file.'
    )
    parser.add_argument(
        '-orfs',
        help='List containing open reading frames as tuples ex:[(2,10), (4,9), (8,0)]',
        default=[(0,100)]
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
