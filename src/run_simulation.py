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
import yaml

from src.sequence_info import NUCLEOTIDES
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
        '--config', default=None,
        help='Path to a YAML file that specifies parameters for the simulation. '
             'Parameters include: mu, kappa, pi, ORF coordinates, and the shape(s) and number(s) '
             'of gamma classes for each ORF as well as the dN and dS values'
    )
    parser.add_argument(
        '--global_rate', type=float, default=1,
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
        '--omega_dist', type=str, default=ss.gamma,
        help='The distribution type from which, from which the omega values are drawn'
    )
    parser.add_argument(
        '--mu_classes', type=int, default=4,
        help='Number of nucleotide classes: The number of classes in which we are going to classify nucleotides on the Event Tree'
    )
    parser.add_argument(
        '--mu_shape', type=float, default=1,
        help='The shape parameter of the mu distribution (log normal or gamma)'
    )
    parser.add_argument(
        '--mu_dist', type=str, default=ss.lognorm,  # ss.lognorm
        help='The distribution type from which get numbers for the mu classes'
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
        orf_list = orf_locations[strand]

        for orf in orf_list:
            orf_length = 0
            orf_start = orf['coords'][0][0]
            orf_end = orf['coords'][-1][1]  # handle spliced ORFs

            orf_length += abs(orf['coords'][-1][1] - orf['coords'][0][0])

            # Check that the start and end positions are integers
            if type(orf['coords'][0][0]) is not int or type(orf['coords'][0][1]) is not int and orf not in invalid_orfs:
                print("Invalid orf: {}; Start and end positions must be integers.".format(orf))
                invalid_orfs[strand].append(orf)

            # Check that the start and stop positions are in the range of the sequence
            if 0 > orf_start or seq_length < orf_start or \
                    0 > orf_end or seq_length < orf_end and orf not in invalid_orfs[strand]:
                print("Invalid orf: {}; Positions must be between 0 and {}".format(orf, seq_length))
                invalid_orfs[strand].append(orf)

            # Check that the ORF range is valid
            if orf_length < 8 and orf not in invalid_orfs[strand]:
                invalid_orfs[strand].append(orf)

            # Check that the ORF is composed of codons
            # Inclusive range (start and end coordinates included)
            if orf_length % 3 != 0 and orf not in invalid_orfs[strand]:
                print("Invalid orf: {}; Not multiple of three".format(orf))
                invalid_orfs[strand].append(orf)

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
            difference = abs(fwd_orf['coords'][0][0] - first_orf['coords'][0][0]) % 3
            if difference == 0:
                sorted_orfs['+0'].append(fwd_orf)
            elif difference == 1:
                sorted_orfs['+1'].append(fwd_orf)
            elif difference == 2:
                sorted_orfs['+2'].append(fwd_orf)

        for rev_orf in reverse_orfs:
            difference = abs(rev_orf['coords'][0][0] - first_orf['coords'][0][0])
            if difference == 0:
                sorted_orfs['-0'].append(rev_orf)
            elif difference == 1:
                sorted_orfs['-1'].append(rev_orf)
            elif difference == 2:
                sorted_orfs['-2'].append(rev_orf)

    return sorted_orfs


def create_values_dict(alpha, ncat, string, dist):
    """
    Creates dictionary with values (as rates) drawn from a discretized gamma distribution
    :param alpha: shape parameter
    :param ncat: Number of categories
    :param dist: the distribution (gamma or log normal)
    :param string: strings like "omega" or "mu" to name the keys of the dictionary
    """

    nt_categories = discretize(alpha, ncat, dist)
    nt_categories_dict = {}
    for i, item in enumerate(nt_categories):
        cat = f"{string}{i+1}"
        nt_categories_dict[cat] = item

    return nt_categories_dict


def discretize(alpha, ncat, dist):
    """
    Divide a distribution into a number of intervals with equal probability and get the mid point of those intervals
    From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
    :param alpha: shape parameter
    :param ncat: Number of categories
    :param dist: distribution of probabilities
    :return: array with ncat number of values
    """

    if dist == ss.gamma:
        dist = dist(alpha, scale=0.4)
    elif dist == ss.lognorm:
        dist = dist(s=alpha, scale=0.5)  # scale=np.exp(0.05 * alpha**2)

    quantiles = dist.ppf(np.arange(0, ncat) / ncat)
    rates = np.zeros(ncat, dtype=np.double)

    for i in range(ncat-1):
        rates[i] = (ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[i], quantiles[i+1])[0])

    rates[ncat-1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x), quantiles[ncat-1], np.inf)[0]

    return rates


def read_sequence(in_file):
    """
    Extract the sequence from the input file
    :param in_file: path to the file containing the sequence.
        Supported file types are GenBank (.gb, .genbank) and FASTA (.fa, .fasta)
    :return: the sequence as a string
    """
    seq = ''
    if in_file.lower().endswith(".gb") or in_file.lower().endswith("genbank"):
        # Loop through records
        for rec in SeqIO.parse(in_file, format="genbank"):
            seq = rec.seq

    elif in_file.lower().endswith(".fasta") or in_file.lower().endswith(".fa"):
        # Read in the sequence
        with open(in_file) as seq_file:
            seq = ''
            for line in seq_file:
                # Skip header if the file is a FASTA file
                if not (line.startswith(">") or line.startswith("#")):
                    seq += line.strip('\n\r').upper()

    else:
        print("Sequence files must end in '.fa', '.fasta', '.gb', 'genbank'")
        logging.error("Invalid file type: files must end in '.fa', '.fasta', '.gb', 'genbank'")
        sys.exit(1)

    return seq


def parse_genbank_orfs(in_seq):
    """
    Extract ORFs from the GenBank file
    """
    orf_locations = {'+': [], '-': []}  # ORF locations sorted by strand

    # Loop through records
    for rec in SeqIO.parse(in_seq, format="genbank"):
        # Read ORFs from GenBank file
        cds = [feat for feat in rec.features if feat.type == "CDS"]
        # Record the first occurrence of the ORFs
        for cd in cds:
            orf = {}
            strand = ''
            for loc in cd.location.parts:
                if loc.strand > 0:
                    strand = '+'
                else:
                    strand = '-'
                orf['coords'] = (int(loc.start), int(loc.end))

            orf_locations[strand].append(orf)

    return orf_locations


def parse_orfs_from_csv(in_orfs, omega_values):
    """
    Reads ORFs from a csv file
    :param in_orfs: orfs specified by the user, default (None)
    :param omega_values: omega values, derived from the discretized gamma distribution
    :return: ORFs as a list of tuples
    """

    orf_locations = {'+': [], '-': []}  # ORF locations sorted by strand
    # Read ORFs from csv file
    with open(in_orfs) as orf_handle:
        for line in orf_handle:
            line = line.strip()
            line = line.split(',')

            for partial_coord in line:
                if ':' in partial_coord:
                    orf = {}
                    partial_coord = partial_coord.split(':')

                    # Check if strand is specified
                    if len(line) == 3:
                        if int(partial_coord[2]) > 0:
                            strand = '+'
                        else:
                            strand = '-'

                    orf['coords'] = (int(partial_coord[0]), int(partial_coord[1]))
                    orf['omega_values'] = omega_values

                    orf_locations[strand].append(orf)
                    continue

                else:
                    orf = {}
                    if int(line[1]) > 0:
                        strand = '+'
                    else:
                        strand = '-'

                    orf['coords'] = (int(line[0]), int(line[1]))
                    orf['omega_values'] = omega_values

            if orf not in orf_locations[strand]:
                orf_locations[strand].append(orf)

    return orf_locations


def parse_orfs_from_yaml(settings):
    """
    Reads ORFs from a YAML file containing the ORF coordinates and the parameters of the dN/dS distribution for each ORF
    :param settings: dictionary representation of YAML file
    :return: orf_locations, a dictionary of ORFs sorted by strand (+ or -)
    """
    orf_locations = {'+': [], '-': []}

    raw_coords = list(settings['orfs'].keys())

    for raw_coord in raw_coords:
        # Spliced ORF
        if ':' in raw_coord:
            orf = {}
            raw_coord = raw_coord.split(':')

            # Read in partial ORFs
            strand = ''
            for coords in raw_coord:
                coords = coords.split(',')
                if len(coords) == 3:
                    if int(coords[2]) > 0:
                        strand = '+'
                    else:
                        strand = '-'
                orf['coords'] = [(int(coords[0]), int(coords[1]))]

                # Get omega values based on the full ORF
                orf['omega_shape'] = settings['orfs'][raw_coord]['omega_shape']
                orf['omega_classes'] = settings['orfs'][raw_coords]['omega_classes']
                orf['omega_values'] = list(discretize(settings['orfs'][raw_coord]['omega_shape'],
                                                      settings['orfs'][raw_coord]['omega_classes'],
                                                      settings['orfs'][raw_coord]['omega_dist']))

                orf_locations[strand].append(orf)

        else:
            orf = {}
            coords = raw_coord.split(',')
            coords = list(map(int, coords))  # Convert string to integer

            orf['coords'] = [coords]
            if coords[0] < coords[1]:
                strand = '+'
            else:
                strand = '-'

            orf['omega_shape'] = settings['orfs'][raw_coord]['omega_shape']
            orf['omega_classes'] = settings['orfs'][raw_coords]['omega_classes']
            orf['omega_values'] = list(discretize(settings['orfs'][raw_coord]['omega_shape'],
                                                  settings['orfs'][raw_coord]['omega_classes'],
                                                  settings['orfs'][raw_coord]['omega_dist']))

            orf_locations[strand].append(orf)

    return orf_locations


def set_global_omega_values(orf_locations, omega_values, omega_shape, omega_classes):
    """
    Sets the dN and dS values for each the reading frames
    :param orf_locations: dictionary of ORFs sorted by the strand
    :param omega_values: list of omega values, derived from the discretized gamma distribution
    return: orf_locations updated to contains dN and dS values for each ORF
    """
    for strand in orf_locations:
        orf_list = orf_locations[strand]
        for orf in orf_list:
            orf['omega_values'] = omega_values
            orf['omega_shape'] = omega_shape
            orf['omega_classes'] = omega_classes

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
    stop = ["TGA", "TAG", "TAA"]
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
        yield my_orf[i:i + 3], i
        i += 3


def count_internal_stop_codons(seq, strand, orf):
    """
    Look for stop codons inside the CDS
    :param seq: the input sequence
    :param strand: the strand (1 or -1)
    :param orf: dictionary containing the coordinates of the ORF and the dN and dS values
    :return: the number of stop codons in the coding sequence
    """
    pat = '(TAA|TGA|TAG)'
    reg = re.compile(pat)
    stop_count, cds = 0, ""

    # Get CDS
    cds += seq[orf['coords'][0]: orf['coords'][1]]

    if strand == '-':    # Reverse strand
        cds = cds[::-1]
    stop_matches = reg.finditer(str(cds))

    # Check if the STOP codon is in the same frame the start codon and is not the last STOP codon in the ORF
    for match in stop_matches:
        stop_start = match.span()[0]
        stop_end = match.span()[1]
        if stop_start % 3 == orf['coords'][0] % 3 and stop_end < orf['coords'][1]:
            stop_count += 1

    return stop_count


def get_pi(pi, settings, s):
    """
    Extracts the value of pi from the command line or a configuration file
    :param pi: the value of pi input from the command line, can be 'None'
    :param settings: a dictionary of settings from the configuration file, can be 'None'
    :param s: the sequence as a string
    :return: the value of pi
    """
    keys = ['A', 'T', 'G', 'C']

    if settings is not None:
        pi = list(settings['pi'].values())
        return dict(zip(keys, pi))

    # If the user specified stationary frequencies
    if all(type(freq) == float for freq in pi):
        pi = dict(zip(keys, pi))
        return pi

    # If the user did not specify stationary frequencies
    elif all(freq is None for freq in pi):
        pi = Sequence.get_frequency_rates(s)
        return pi

    else:
        print("Invalid input: {}".format(pi))
        sys.exit(1)


def get_kappa(kappa, settings):
    if settings:
        kappa = settings['kappa']
    return kappa


def get_global_rate(global_rate, settings):
    if settings:
        global_rate = settings['global_rate']
    return global_rate


def get_mu_classes(mu_classes, settings):
    if settings:
        mu_classes = settings['mu_classes']
    return mu_classes


def get_mu_shape(mu_shape, settings):
    if settings:
        mu_shape = settings['mu_shape']
    return mu_shape


def get_mu_dist(mu_dist, settings):
    if settings:
        mu_dist = settings['mu_dist']
    return mu_dist


def read_settings_from_config(config_path):
    """
    Parse the config file if it exists
    :param config_path: path to the configuration file
    :return: a dictionary of settings from the configuration file (if no config file exists, the dictionary is empty)
    """
    settings = {}

    if config_path:
        with open(config_path, 'r') as stream:
            try:
                settings = yaml.safe_load(stream)
            except yaml.YAMLError as e:
                print(e)

    return settings


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
    LOG_FILENAME = f'{file_name.split(".")[0]}_evol_simulation.log'
    logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
    logging.info(f"\nSimulation started at: {start_time}")

    # Read in sequence
    s = read_sequence(args.seq)

    # Parse configuration file (if it exists)
    if args.config:
        with open(args.config, 'r') as stream:
            try:
                settings = yaml.safe_load(stream)
            except yaml.YAMLError as e:
                print(e)
    else:
        settings = {}

    # Get parameters for run
    pi = get_pi(args.pi, settings, s)
    global_rate = get_global_rate(args.global_rate, settings)
    kappa = get_kappa(args.kappa, settings)
    mu_classes = get_mu_classes(args.mu_classes, settings)
    mu_shape = get_mu_shape(args.mu_shape, settings)
    mu_dist = get_mu_dist(args.mu_dist, settings)

    # Create classes to classify nucleotides on the Event Tree
    mu_values = create_values_dict(mu_shape, mu_classes, "mu", mu_dist)
    print(f"Categories: {mu_values}")

    # Log global parameters
    logging.info(f"Parameters for the run: \n"
                 f"Pi: {pi}\n"
                 f"Global rate: {global_rate}\n"
                 f"Kappa: {kappa}\n"
                 f"Number of nucleotide classification classes: {mu_classes}\n"
                 f"Nucleotide classification shape parameter: {mu_shape}\n"
                 f"Rates classification values: {mu_values}")

    # Use a different omega distribution for each ORF
    if args.config:
        orf_locations = parse_orfs_from_yaml(settings)

    # Use the same omega distribution for all ORFs
    else:
        omega_values = list(discretize(args.omega_shape, args.omega_classes, args.omega_dist))

        # Read ORFs from GenBank file
        if args.seq.lower.endswith('.gb') or args.seq.lower.endswith('genbank'):
            orfs = parse_genbank_orfs(args.seq)

        # Read in ORFs from .csv file
        if args.orfs:
            orfs = parse_orfs_from_csv(args.orfs, omega_values)

        orf_locations = set_global_omega_values(orfs, omega_values, args.omega_shape, args.omega_classes)

    # Check if the ORFs are valid
    invalid_orfs = valid_orfs(orf_locations, len(s))

    # Omit the invalid ORFs
    if invalid_orfs['+'] or invalid_orfs['+']:
        invalid_orf_msg = ""
        for strand in invalid_orfs:
            orfs = invalid_orfs[strand]
            for orf in orfs:
                invalid_orf_msg += f" {orf} "
                orf_locations[strand].remove(orf)

        if invalid_orf_msg:
            print(f"\nOmitted orfs: {invalid_orf_msg}\n")
            logging.warning(f"Omitted orfs: {invalid_orf_msg}")

    # Check if the CDSs have stop codons inside them
    for strand in orf_locations:
        orfs = orf_locations[strand]
        for orf_coords in orfs:
            if strand == '-':       # Reverse strand
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

    # Record the parameters for each ORF
    for frame in orfs:
        orf_list = orfs[frame]
        for orf in orf_list:
            for c in orf['coords']:
                coords = ','.join(map(str, c))

    # Read in the tree
    phylo_tree = Phylo.read(args.tree, 'newick', rooted=True)
    # logging.info("Phylogenetic tree: {}".format(args.tree))

    # # Make Sequence object
    print("\nCreating root sequence")
    root_sequence = Sequence(s, orfs, kappa, global_rate, pi, mu_values, args.circular)

    # Run simulation
    # print("\nRunning simulation")
    simulation = SimulateOnTree(root_sequence, phylo_tree, args.outfile)
    simulation.get_alignment(args.outfile)

    end_time = datetime.now()
    print(f"Simulation duration: {end_time - start_time} seconds")
    logging.info(f"\nSimulation Ended at: {end_time}\n"
                 f"Simulation lasted: {end_time - start_time} seconds\n")


if __name__ == '__main__':
    main()
