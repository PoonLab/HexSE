import argparse
import re
import sys
import logging

import numpy as np
import scipy
import scipy.stats as ss
from Bio import Phylo
from Bio import SeqIO
from datetime import datetime

from src.sequence_info import NUCLEOTIDES, COMPLEMENT_DICT
from src.sequence_info import Sequence
from src.simulation import SimulateOnTree


def get_args(parser):
    parser.add_argument(
        'seq',
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'tree',
        help='Path to file containing phylogenetic tree in Newick format.'
    )
    parser.add_argument(
        '--outfile', default=None, help='Path to the alignment file.'
    )
    parser.add_argument(
        '--orfs', default=None,
        help='Path to a csv file containing the start and end coordinates of the open reading frames. '
             'Format: start,end'
             'If no ORFS are specified, the program will find ORFs automatically'
    )
    parser.add_argument(
        '--mu', type=float, default=0.0005,
        help='Global substitution rate per site per unit time'
    )
    parser.add_argument(
        '--kappa', type=float, default=0.3,
        help='Transversion/ transition rate assuming time reversibility.'
    )
    parser.add_argument(
        '--pi', type=float, default=[None, None, None, None],
        help='Vector of stationary nucleotide frequencies. If no value is specified, '
             'the program will use the empirical frequencies in the sequence. Format: [A, T, G, C]'
    )
    parser.add_argument(
        '--omega_classes', type=int, default=4,
        help='The number of omega classes'
    )
    parser.add_argument(
        '--omega_shape', type=float, default=2.0,
        help='The shape parameter of the gamma distribution, from which omega values are drawn'
    )
    parser.add_argument(
        '--nc_classes', type=int, default=4,
        help='Number of nucleotide classes: The number of classes in which we are going to classify nucleotides on the Event Tree'
    )
    parser.add_argument(
        '--nc_shape', type=float, default=0.1,
        help='The shape parameter of the gamma distribution, from which rate values are drawn to clasify nucleotides on the Event Tree'
    )
    parser.add_argument(
        '--circular', action='store_true',
        help='True for circular genomes. By default, false for linear genomes'
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


def valid_orfs(orf_locations, seq_length):
    """
    Verifies that the input ORFs are a list of tuples containing the start and end positions of ORFs.
    Example of valid input: [(1, 9), (27, 13)]
    Example of invalid input: (1, 9), (27, 13)
    :param orf_locations: The list of open reading frames
    :param seq_length: The length of the original sequence
    :return: <True> if the ORFs are valid, <False> otherwise
    """
    invalid_orfs = {'+': [], '-': []}   # ORF locations sorted by strand

    for strand in orf_locations:
        orf_coords = orf_locations[strand]

        for orf_coord in orf_coords:
            orf_length = 0
            orf_start = orf_coord[0][0]
            orf_end = orf_coord[-1][1]

            for coord in orf_coord:
                orf_length += abs(coord[1] - coord[0])

                # Check that the start and end positions are integers
                if type(coord[0]) is not int or type(coord[1]) is not int and orf_coord not in invalid_orfs:
                    print("Invalid orf: {}; Start and end positions must be integers.".format(orf_coord))
                    invalid_orfs[strand].append(orf_coord)

            # Check that the start and stop positions are in the range of the sequence
            if 0 > orf_start or seq_length < orf_start or \
                    0 > orf_end or seq_length < orf_end and orf_coord not in invalid_orfs[strand]:
                print("Invalid orf: {}; Positions must be between 0 and {}".format(orf_coord, seq_length))
                invalid_orfs[strand].append(orf_coord)

            # Check that the ORF range is valid
            if orf_length < 8 and orf_coord not in invalid_orfs[strand]:
                invalid_orfs[strand].append(orf_coord)

            # Check that the ORF is composed of codons
            # Inclusive range (start and end coordinates included)
            if orf_length % 3 != 0 and orf_coord not in invalid_orfs[strand]:
                print("Invalid orf: {}; Not multiple of three".format(orf_coord))
                invalid_orfs[strand].append(orf_coord)

    return invalid_orfs


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
                    orf = (position, match.end())
                    reading_frames.append(orf)
                    break

    # Forward (positive) reading frames of the reverse complement of the original
    # sequence is equivalent to reverse (negative) reading frames of the original sequence
    rcseq = Sequence.complement(seq, rev=True)

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
                    orf = (len(rcseq) - position, len(rcseq) - match.end())
                    reading_frames.append(orf)
                    break

    return reading_frames


def sort_orfs(orf_locations):
    """
    Store ORFs in position according to plus zero ORF (first of the list).
    They will be classified as (+0, +1, +2, -0, -1, -2)
    :return sorted_orfs: List of ORFs classified according to their shift relative to the
                        plus zero reading frame (+0, +1, +2, -0, -1, -2)
    """
    sorted_orfs = {'+0': [], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}

    if orf_locations['+']:

        forward_orfs = orf_locations['+']
        reverse_orfs = orf_locations['-']
        first_orf = forward_orfs[0]

        for fwd_orf in forward_orfs:
            difference = abs(fwd_orf[0][0] - first_orf[0][0]) % 3
            if difference == 0:
                sorted_orfs['+0'].append(fwd_orf)
            elif difference == 1:
                sorted_orfs['+1'].append(fwd_orf)
            elif difference == 2:
                sorted_orfs['+2'].append(fwd_orf)

        for rev_orf in reverse_orfs:
            difference = abs(rev_orf[0][0] - first_orf[0][0])
            if difference == 0:
                sorted_orfs['-0'].append(rev_orf)
            elif difference == 1:
                sorted_orfs['-1'].append(rev_orf)
            elif difference == 2:
                sorted_orfs['-2'].append(rev_orf)

    return sorted_orfs


def get_rate_values(alpha, ncat):
    """
    Draw ncat number of dN or dS values from a discretized gamma distribution
    :param alpha: shape parameter
    :param ncat: Number of categories (expected dS values)
    :return: list of ncat number of omega values (e.i. if ncat = 3, omega_values = [0.29, 0.65, 1.06])

    """
    values = discretize_gamma(alpha, ncat)
    rate_values = list(values)
    return rate_values


def discretize_gamma(alpha, ncat):
    """
    Divide the gamma distribution into a number of intervals with equal probability and get the mid point of those intervals
    From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
    :param alpha: shape parameter
    :param ncat: Number of categories
    :return: array with ncat number of values
    """
    dist = ss.gamma(alpha, scale=1 / alpha)

    quantiles = dist.ppf(np.arange(0, ncat) / ncat)
    rates = np.zeros(ncat, dtype=np.double)  # return a new array of shape ncat and type double

    for i in range(ncat - 1):
        rates[i] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[i], quantiles[i + 1])[0]

    return rates


def parse_genbank(in_seq, in_orfs=None):
    """
    When input is in <genbank> format, extract nucleotide sequence and orfs (in case user does not specify).
    :param in_seq: sequence in genbank Format
    :param in_orfs: file handle containing the orfs
    :return tuple: (sequence, orfs)
    """
    orf_locations = {'+': [], '-': []}  # ORF locations sorted by strand
    # Loop through records
    for rec in SeqIO.parse(in_seq, format="genbank"):
        seq = rec.seq

        if in_orfs is None:  # User did not specify ORFs
            cds = [feat for feat in rec.features if feat.type == "CDS"]
            # Record the first occurrence of the ORFs
            for cd in cds:
                coords = []
                strand = ''
                for loc in cd.location.parts:
                    if loc.strand > 0:
                        strand = '+'
                    else:
                        strand = '-'
                    coords.append((int(loc.start), int(loc.end)))
                orf_locations[strand].append(coords)

        else:
            orf_locations = check_orfs(in_orfs)

    return seq, orf_locations


def parse_fasta(in_seq):
    """
    If input is a fasta file, retrieve nucleotide sequence
    :param in_seq: the sequence
    :return s: the nucleotide sequence
    """
    # Read in the sequence
    with open(in_seq) as seq_file:
        s = ''
        for line in seq_file:
            # Skip header if the file is a FASTA file
            if not (line.startswith(">") or line.startswith("#")):
                s += line.strip('\n\r').upper()
    return s


def check_orfs(in_orfs=None, s=None):
    """
    :param in_orfs: orfs specified by the user, default (None)
    :param s: the original sequence
    :return: ORFs as a list of tuples
    """
    # Check if the user specified orfs
    if in_orfs is None:
        orf_locations = get_open_reading_frames(s)

    # Read ORFs as a list of tuples
    else:
        orf_locations = {'+': [], '-': []}  # ORF locations sorted by strand
        with open(in_orfs) as orf_handle:
            for line in orf_handle:
                line = line.strip()

                # Spliced ORF
                if ':' in line:
                    orf_coords = []
                    line = line.split(':')

                    # Read in partial ORFs
                    strand = ''
                    for coords in line:
                        coords = coords.split(',')
                        if len(coords) == 3:
                            if int(coords[2]) > 0:
                                strand = '+'
                            else:
                                strand = '-'
                        orf = (int(coords[0]), int(coords[1]))
                        orf_coords.append(orf)

                    # Check if the strand is valid
                    if strand != '+' and strand != '-':
                        print("Invalid strand: {}".format(strand))
                        sys.exit(1)

                    orf_locations[strand].append(orf_coords)

                else:
                    line = line.split(',')
                    if int(line[2]) > 0:
                        strand = '+'
                    else:
                        strand = '-'

                    # Check if the strand is valid
                    if strand != '+' and strand != '-':
                        print("Invalid strand: {}".format(strand))
                        sys.exit(1)

                    orf = (int(line[0]), int(line[1]))
                    orf_locations[strand].append([orf])

    return orf_locations


def create_log_file(input_file_name):
    """
    Create a log file with information for the run
    """
    # Select name of the input file
    file_name = input_file_name.split("/")[-1]
    file_name = file_name.split(".")[0]
    LOG_FILENAME = "{}_evol_simulation.log".format(file_name)

    return LOG_FILENAME


def stop_in_seq(seq, start, end):
    """
    Look for stop codons inside the CDS
    """
    cds = seq[start:end]
    stop =  ["TGA", "TAG", "TAA"]
    stop_count = 0
    for codon, nt in codon_iterator(cds, start, end):
        if codon in stop:
            stop_count += 1

    return stop_count


def codon_iterator(my_orf, start_pos, end_pos):
    """
    Generator to move every three nucleotides (codon)
    :param my_orf: A list of Nucleotides in the ORF
    :param start_pos: The start position of the ORF
    :param end_pos: The end position of the ORF
    :yield codon
    """
    if start_pos > end_pos:  # Negative strand
        my_orf.reverse()
    i = 0
    while i < len(my_orf):
        yield (my_orf[i:i + 3], i)
        i += 3


def count_internal_stop_codons(seq, strand, orf_coords):
    """
    Look for stop codons inside the CDS
    :param seq: the input sequence
    :param strand: the strand (1 or -1)
    :param orf_coords: list containing the coordinates of the ORF and the frame
    :return: the number of stop codons in the coding sequence
    """
    pat = '(TAA|TGA|TAG)'
    reg = re.compile(pat)
    stop_count, cds = 0, ""

    # Get CDS
    for coord in orf_coords:
        cds += seq[coord[0]:coord[1]]

    if strand == '-':    # Reverse strand
        cds = cds[::-1]
    stop_matches = reg.finditer(str(cds))

    # Check if the STOP codon is in the same frame the start codon and is not the last STOP codon in the ORF
    for match in stop_matches:
        stop_end = match.span()[1]
        if (stop_end - 1) % 3 == orf_coords[0][0] % 3 and stop_end > orf_coords[-1][1]:
            stop_count += 1

    return stop_count

def create_nucleotide_categories_dict(alpha, ncat):
    """
    Creates categories drawn from a discretized gamma distribution to store each nucleotide according to its mutation rates on the Event Tree
    :param alpha: shape parameter
    :param ncat: Number of catefories
    """

    nt_categories = get_rate_values(alpha, ncat)
    nt_categories_dict = {}
    for i, item in enumerate(nt_categories):
        cat = f"cat{i+1}"
        nt_categories_dict[cat]={'value':item}

    return nt_categories_dict

def main():
    start_time = datetime.now()
    print("\nStarted at: ", datetime.now())

    parser = argparse.ArgumentParser(
        description='Simulates and visualizes the evolution of a sequence through a phylogeny'
    )
    args = get_args(parser)
    input = args.seq.lower()

    # Create log file
    file_name = input.split("/")[-1]
    LOG_FILENAME = "{}_evol_simulation.log".format(file_name.split(".")[0])
    logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)

    # Check input format for the nucleotide sequence
    if input.endswith(".gb") or input.endswith("genbank"):  # If genbank file
        s, orf_locations = parse_genbank(args.seq, args.orfs)
        logging.info("Input sequence is: {} in genbank format".format(input))

    elif input.endswith(".fasta") or input.endswith(".fa"):   # If fasta file
        s = parse_fasta(args.seq)
        orf_locations = check_orfs(args.orfs, s)
        logging.info("Input sequence is: {} in fasta format".format(input))

    else:
        print("Sequence files must end in '.fa', '.fasta', '.gb', 'genbank'")
        logging.error("Invalid sequence: files must end in '.fa', '.fasta', '.gb', 'genbank'")
        sys.exit()

    # Check if the ORFs are valid
    invalid_orfs = valid_orfs(orf_locations, len(s))

    # Omit the invalid ORFs
    if invalid_orfs:
        invalid_orf_msg = ""
        for strand in invalid_orfs:
            orfs = invalid_orfs[strand]
            for orf in orfs:
                invalid_orf_msg += " {} ".format(orf)
                orf_locations[strand].remove(orf)

        if invalid_orf_msg:
            print("\nOmitted orfs: {}\n".format(invalid_orf_msg))
            logging.warning("Omitted orfs: {}".format(invalid_orf_msg))

    # Check if the CDSs have stop codons inside them
    for strand in orf_locations:
        orfs = orf_locations[strand]
        for orf_coords in orfs:
            if strand == -1:       # Reverse strand
                stop_count = count_internal_stop_codons(Sequence.complement(s), strand, orf_coords)
            else:                   # Forward strand
                stop_count = count_internal_stop_codons(s, strand, orf_coords)

            # CDS has more than one stop codon (the final one)
            if stop_count > 1:
                orf_locations[strand].remove(orf_coords)
                print(f"Omitted orf: {orf_coords} has {stop_count} STOP codons")

    # Since ORFs are valid, sort the ORFs by reading frame
    orfs = sort_orfs(orf_locations)
    logging.info("Valid orfs: {}".format(orfs))

    # Check if sequence is valid
    if not valid_sequence(s):
        print("Invalid sequence: {}".format(s))
        logging.error("Invalid sequence: {}".format(s))
        sys.exit(0)

    # If the user did not specify stationary frequencies
    if all(freq is None for freq in args.pi):
        pi = Sequence.get_frequency_rates(s)

    # If the user specified stationary frequencies
    elif all(freq is type(float) for freq in args.pi):
        keys = ['A', 'T', 'G', 'C']
        pi = dict(zip(keys, args.pi))

    else:
        print("Invalid input: {}".format(args.pi))
        exit(0)

    # Draw omeg values and create classes to classify nucleotides on the Event Tree
    omega_values = get_rate_values(args.omega_shape, args.omega_classes)
    nt_categories_dict = create_nucleotide_categories_dict(args.nc_shape, args.nc_classes)

    logging.info(f"Parameters for the run: \nPi: {pi}\nMu: {args.mu}\nKappa: {args.kappa}\nNumber of omega classes: {args.omega_classes}\n\
    Omega shape parameter: {args.omega_shape}\nRates classification values: {nt_categories_dict}\n\
    Number of nucleotide classification classes: {args.nc_classes}\nNucleotide clasification shape parameter: {args.nc_shape}")

    # TODO: Allow user to modify alpha and ncat. Use this values as default

    print(">>> CATEGORIES")
    print(nt_categories_dict)

    # Read in the tree
    # phylo_tree = Phylo.read(args.tree, 'newick', rooted=True)
    # logging.info("Phylogenetic tree: {}".format(args.tree))
    #
    # # Make Sequence object
    print("\nCreating root sequence")
    root_sequence = Sequence(s, orfs, args.kappa, args.mu, pi, omega_values, nt_categories_dict, args.circular)


    # # Run simulation
    # print("\nRunning simulation")
    # simulation = SimulateOnTree(root_sequence, phylo_tree, args.outfile)
    # simulation.get_alignment(args.outfile)
    #
    # print("Simulation duration: {} seconds".format(datetime.now() - start_time))


if __name__ == '__main__':
    main()
