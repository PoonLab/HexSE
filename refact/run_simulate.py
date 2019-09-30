import argparse
import sys
import re

from Bio import Phylo

from refact.new_sequence_info import NUCLEOTIDES
from refact.new_sequence_info import COMPLEMENT_DICT
from refact.new_sequence_info import Sequence


def get_args(parser):
    parser.add_argument(
        'seq', type=argparse.FileType('r'),
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'tree',
        help='Path to file containing phylogenetic tree in Newick format.'
    )
    parser.add_argument(
        'outfile', type=argparse.FileType('w'),
        help='Path to the alignment file.'
    )
    parser.add_argument(
        '--orfs', type=argparse.FileType('r'), default=None,
        help='Path to a csv file containing the start and end coordinates of the open reading frames'
    )
    parser.add_argument(
        '--mu', type=float,
        help='Global substitution rate per site per unit time'
    )
    parser.add_argument(
        '--kappa', nargs="6", type=float, default=[1, 1, 1, 1, 1, 1],
        help='List of transversion/ transition rates assuming time reversibility. '
             'Format: [AC, AG, AT, CG, CT, GT]'
    )
    parser.add_argument(
        '--pi', nargs="4", type=float, default=[None, None, None, None],
        help='Vector of stationary nucleotide frequencies. If no value is specified, '
             'the program will use the empirical frequencies in the sequence. Format: [A, T, G, C]'
    )
    parser.add_argument(
        '--omega', nargs="6", type=float, default=[None, None, None, None, None, None],
        help='List of dN/dS ratios for each reading frame along the length of the sequence. '
             'Format: [+0, +1, +2, -0, -1, -2]'
    )

    return parser.parse_args()


def valid_sequence(seq):
    """
     Verifies that the length of the input sequence is valid and the sequence is composed of only nucleotides.
    A valid sequence is assumed to be composed of a START codon, at least one amino acid codon, and a STOP codon.
    :return is_valid: <True> if the sequence is valid, <False> otherwise
    """
    is_valid = len(seq) >= 9 and all(pos in NUCLEOTIDES for pos in seq.upper())
    return is_valid


def valid_orfs(orfs, seq):
    """
    Verifies that the input ORFs are a list of tuples containing the start and end positions of ORFs.
    Example of valid input: [(1, 9), (27, 13)]
    Example of invalid input: (1, 9), (27, 13)
    :param orfs: The list of open reading frames
    :param seq: The original sequence as a string
    :return: <True> if the ORFs are valid, <False> otherwise
    """
    if not type(orfs) == list:
        print("Invalid format: {} \nOpen reading frames must of the following format: [(0, 8), ...] where 0  "
              "is the start position of the orf, and 8 is the end position.".format(orfs))
        return False

    for orf in orfs:
        # Check that each orfs is a tuple of length 2
        if type(orf) is not tuple or len(orf) is not 2:
            print("Invalid orf: {} \nOpen reading frames must of the following format: [(0, 8), ...] where 0  "
                  "is the start position of the orf, and 8 is the end position.".format(orf))
            return False

        # Check that the ORF range is valid
        if orf[0] == orf[1]:
            print("Invalid orf: {}".format(orf))
            return False

        # Check that the start and end positions are integers
        if type(orf[0]) is not int or type(orf[1]) is not int:
            print("Invalid orf: {} \nStart and end positions must be integers.".format(orf))
            return False

        # Check that the start and stop positions are in the range of the sequence
        if 0 > orf[0] or len(seq) < orf[0] or 0 > orf[1] or len(seq) < orf[1]:
            print("Invalid orf: {} \nPositions must be between 0 and {}".format(orf, len(seq)))
            return False

        # Check that the ORF is composed of codons
        if orf[1] > orf[0]:  # Forward strand
            if (orf[1] - orf[0]) % 3 != 2:
                print("Invalid orf: {}\n ORFs must be composed of codons".format(orf))
                return False
        if orf[0] > orf[1]:  # Reverse strand
            if (orf[0] - orf[1]) % 3 != 2:
                print("Invalid orf: {}\n ORFs must be composed of codons".format(orf))
                return False

    return True


def reverse_and_complement(seq):
    """
    Generates the reverse complement of a DNA sequence
    :param: my_region <option> A sub-sequence of the original sequence
    :return rcseq: The reverse complement of the sequence
    """
    rseq = reversed(seq.upper())
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += COMPLEMENT_DICT[i]
    return rcseq


def get_open_reading_frames(seq):
    """
    Gets positions of the START and STOP codons for each open reading frame in the forward and reverse directions.
    Positions of the START and STOP codons indexed relative to the forward strand.
    :return reading_frames: a list of tuples containing the index of the first nucleotide of
                the START codon and the index of the last nucleotide of the STOP codon
    """
    start_codon = re.compile('ATG', flags=re.IGNORECASE)
    stop = re.compile('(TAG)|(TAA)|(TGA)', flags=re.IGNORECASE)
    reading_frames = []

    # Record positions of all potential START codons in the forward (positive) reading frame
    fwd_start_positions = [match.start() for match in start_codon.finditer(seq)]

    # Find open forward open reading frames
    for position in fwd_start_positions:
        frame = position % 3

        internal_met = False
        # If the ATG codon is an internal methionine and not an initiation codon
        for orf in reversed(reading_frames):

            # If the START codon and the potential START codon are in the same reading frame
            # and the existing ORF ends before the potential ORF, stop searching
            if orf[0] % 3 == frame and orf[1] < position:
                break

            # If the potential START codon is between the range of the START and STOP codons,
            # and it is in the same frame, the codon is an internal methionine
            if orf[0] < position < orf[1] and orf[0] % 3 == frame:
                internal_met = True
                break

        # If the ATG is a START codon and not simply methionine
        if not internal_met:
            for match in stop.finditer(seq, position):
                orf_length = match.end() - position
                # Find a stop codon and ensure ORF length is sufficient in the forward strand
                if match.start() % 3 == frame and orf_length >= 8:
                    # Get the positions in the sequence for the first and last nt of the RF
                    orf = (position, match.end() - 1)
                    reading_frames.append(orf)
                    break

    # Forward (positive) reading frames of the reverse complement of the original
    # sequence is equivalent to reverse (negative) reading frames of the original sequence
    rcseq = reverse_and_complement(seq)

    # Record positions of all potential START codons in the reverse (negative) reading frame
    rev_start_positions = [match.start() for match in start_codon.finditer(rcseq)]

    # Find reverse open reading frames
    for position in rev_start_positions:
        frame = position % 3

        internal_met = False
        # If the ATG codon is an internal methionine and not an initiation codon
        for orf in reversed(reading_frames):

            # If the START codon and the potential START codon are in the same reading frame
            # and the existing ORF ends before the potential ORF, stop searching
            if orf[0] % 3 == frame and orf[1] < position:
                break

            # If the potential START codon is between the range of the START and STOP codons,
            # and it is in the same frame, the codon is an internal methionine
            if orf[0] < position < orf[1] and orf[0] % 3 == frame:
                internal_met = True
                break

        # If the ATG is a START codon and not simply methionine
        if not internal_met:
            for match in stop.finditer(rcseq, position):
                orf_length = match.end() - position
                # Find a stop codon and ensure ORF length is sufficient in the forward strand
                if match.start() % 3 == frame and orf_length >= 8:
                    # Get the positions in the sequence for the first and last nt of the RF
                    orf = (len(rcseq) - 1 - position, len(rcseq) - match.end())
                    reading_frames.append(orf)
                    break

    return reading_frames


def sort_orfs(unsorted_orfs):
    """
    Store ORFs in position according to plus zero ORF (first of the list).
    They will be classified as (+0, +1, +2, -0, -1, -2)
    :return sorted_orfs: List of ORFs classified according to their shift relative to the
                        plus zero reading frame (+0, +1, +2, -0, -1, -2)
    """
    sorted_orfs = {'+0': [], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}

    if unsorted_orfs:
        first_orf = unsorted_orfs[0]
        for orf in unsorted_orfs:
            difference = abs(orf[0] - first_orf[0]) % 3

            if first_orf[0] < first_orf[1]:
                if orf[0] < orf[1]:  # positive strand
                    if difference == 0:
                        sorted_orfs['+0'].append(orf)
                    elif difference == 1:
                        sorted_orfs['+1'].append(orf)
                    elif difference == 2:
                        sorted_orfs['+2'].append(orf)

                elif orf[0] > orf[1]:  # negative strand
                    if difference == 0:
                        sorted_orfs['-2'].append(orf)
                    elif difference == 1:
                        sorted_orfs['-1'].append(orf)
                    elif difference == 2:
                        sorted_orfs['-0'].append(orf)

            else:
                if orf[0] < orf[1]:  # positive strand
                    if difference == 0:
                        sorted_orfs['+2'].append(orf)
                    elif difference == 1:  # plus one
                        sorted_orfs['+1'].append(orf)
                    elif difference == 2:  # plus two
                        sorted_orfs['+0'].append(orf)

                elif orf[0] > orf[1]:  # negative strand
                    if difference == 0:
                        sorted_orfs['-0'].append(orf)
                    elif difference == 1:
                        sorted_orfs['-1'].append(orf)
                    elif difference == 2:
                        sorted_orfs['-2'].append(orf)

    return sorted_orfs


def main():
    parser = argparse.ArgumentParser(
        description='Simulates and visualizes the evolution of a sequence through a phylogeny'
    )
    args = get_args(parser)

    # Read in the sequence
    s = ''
    for line in args.seq:
        # Skip header if the file is a FASTA file
        if not (line.startswith(">") or line().startswith("#")):
            s += line.strip('\n\r').upper()

    # Check if sequence is valid
    if not valid_sequence(s):
        print("Invalid sequence: {}".format(s))
        sys.exit(0)

    # Read in ORFs as a list of tuples
    unsorted_orfs = []
    for line in args.orfs:
        line = line.split(',')
        orf = (line[0], line[1])
        unsorted_orfs.append(orf)

    # Check if ORFs are valid
    if unsorted_orfs is not None:
        if not valid_orfs(unsorted_orfs, s):
            sys.exit(0)
        # Since ORFs are valid, sort the ORFs by reading frame
        orfs = sort_orfs(unsorted_orfs)

    # If the user did not specify ORFs
    else:
        unsorted_orfs = get_open_reading_frames(s)
        orfs = sort_orfs(unsorted_orfs)

    # Read in the tree
    tree = Phylo.read(args.tree, 'newick', rooted=True)

    # Check that all stationary frequencies are non-None
    if args.pi is not None:
        for value in args.pi:
            if value is None:
                print("Invalid input: Missing stationary frequencies.")
                sys.exit(0)

            # Since input pi is valid, reformat the values into a dictionary
            keys = ['A', 'T', 'G', 'C']
            pi = dict(zip(keys, args.pi))

    # If the user does not specify pi values
    else:
        pi = Sequence.get_frequency_rates(s)

    # Make Sequence object
    seq = Sequence(s, reverse_and_complement(s), orfs, args.mu, pi, args.kappa)

    # Reformat omega into a dictionary
    # TODO: check that there is a check for all non-None values in omega
    # keys = ['+0', '+1', '+2', '-0', '-1', '-2']
    # omega = dict(zip(keys, args.omega))

    # Run simulation
    # sim = Simulate(rates, tree)
    #
    # sim.traverse_tree()
    # sim.get_alignment(args.outfile)


if __name__ == '__main__':
    main()
