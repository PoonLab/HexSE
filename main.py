from ovrf_functions import get_reading_frames
from ovrf_functions import sort_orfs
from evol_rates import NUCLEOTIDES
from sequence_info import Sequence
from sequence_info import Nucleotide

import argparse


def valid_sequence(seq):
    """
    Verifies that the length of the input sequence is valid and the sequence is composed of only nucleotides.
    Note: Assumes a valid sequence is composed of a START codon, at least one amino acid codon, and a STOP codon.
    :param seq: the input sequence
    :return valid: true if the sequence is valid, false otherwise
    """
    is_valid = len(seq) >= 9 and all(pos in NUCLEOTIDES for pos in seq)

    if not is_valid:
        raise ValueError("Invalid sequence: {}".format(seq))

    return is_valid


def get_args(parser):
    parser.add_argument(
        'infile',
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'outfile',
        help='Path to the output file'
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

    if valid_sequence(s):
        orf_positions = get_reading_frames(s)

        if orf_positions:
            orfs = sort_orfs(orf_positions)


if __name__ == '__main__':
    main()
