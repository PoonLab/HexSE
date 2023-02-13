import argparse
import re
import sys
import logging
import pprint
import random
import json

import numpy as np
import scipy
import scipy.stats as ss
from Bio import Phylo
from Bio import SeqIO
from datetime import datetime

from .sequence_info import NUCLEOTIDES
from .sequence_info import AMBIGUOUS_NUCLEOTIDES
from .sequence_info import Sequence
from .simulation import SimulateOnTree
from .settings import Settings
from .discretize import discretize


def get_args(parser: argparse.ArgumentParser) -> argparse.Namespace:
    # positional arguments (required)
    parser.add_argument(
        'seq',
        help='Path to the file containing the query sequence.'
    )
    parser.add_argument(
        'tree',
        help='Path to file containing phylogenetic tree in Newick format.'
    )
    parser.add_argument(
        'config',
        help='Path to a YAML file that specifies parameters for the simulation. '
             'Parameters include: mu, kappa, pi, ORF coordinates, and the shape(s) and number(s) '
             'of gamma classes for each ORF as well as the dN and dS values'
    )

    # keyword arguments
    parser.add_argument(
        '--outfile', default=None, help='Path to the alignment file; defaults to stdout.'
    )

    parser.add_argument(
        '--logfile', default=None, help='Path to the log file; defaults to stdout.'
    )

    parser.add_argument(
        '--ci', action='store_true', help='optional, return codon info.'
    )

    return parser.parse_args()


def valid_sequence(seq: str) -> bool:
    """
    Verifies that the length of the input sequence is valid and the sequence is composed of only nucleotides.
    A valid sequence is assumed to be composed of a START codon, at least one amino acid codon, and a STOP codon.
    :return is_valid: <True> if the sequence is valid, <False> otherwise
    """
    is_valid = len(seq) >= 9 and all(pos in NUCLEOTIDES for pos in seq.upper())
    return is_valid


def resolve_ambiguities(seq: str) -> str:
    """
    Resolves ambiguous positions in nt seq by randomly selecting one of its posibilites
    :return: seq without ambiguities
    """
    new_seq = list(seq.upper())
    for id_nt, nt in enumerate(new_seq):
        unambiguous_list = AMBIGUOUS_NUCLEOTIDES.get(nt)
        if unambiguous_list:
            new_seq[id_nt] = random.choice(unambiguous_list)

    return("".join(new_seq))


def valid_orfs(orf_locations: dict, seq_length: int) -> (dict, list):
    """
    Verifies that the input ORFs are a list of tuples containing the start and end positions of ORFs.
    Example of valid input: [(1, 9), (27, 13)]
    Example of invalid input: (1, 9), (27, 13)
    :param orf_locations: The list of open reading frames
    :param seq_length: The length of the original sequence
    :return: <True> if the ORFs are valid, <False> otherwise
    """
    invalid_orfs = {'+': [], '-': []}   # ORF locations sorted by strand
    invalid_msg = []  # Store the reason why and orf was invalid

    for strand in orf_locations:
        orf_list = orf_locations[strand]  # List of dictionaries with ORF info

        for orf in orf_list:
            
            orf_length = 0
            orf_start = orf['coords'][0][0]
            orf_end = orf['coords'][-1][1]  # handle spliced ORFs
            orf_length += abs(sum([end - start for end, start in orf['coords']]))  # Get orf length by adding all fragments

            # Check that the start and end positions are integers
            if type(orf['coords'][0][0]) is not int or type(orf['coords'][0][1]) is not int and orf not in invalid_orfs:
                print("Invalid orf: {}; Start and end positions must be integers.".format(orf))
                invalid_orfs[strand].append(orf)
                invalid_msg.append([orf['coords'],"Start and end position must be intergers"])

            # Check that the start and stop positions are in the range of the sequence
            if 0 > orf_start or seq_length < orf_start or \
                    0 > orf_end or seq_length < orf_end and orf not in invalid_orfs[strand]:
                print("Invalid orf: {}; Positions must be between 0 and {}".format(orf, seq_length))
                invalid_orfs[strand].append(orf)
                invalid_msg.append([orf['coords'],"Positions must be between 0 and sequence length"])

            # Check that the ORF range is valid
            if orf_length < 8 and orf not in invalid_orfs[strand]:
                invalid_orfs[strand].append(orf)

            # Check that the ORF is composed of codons
            # Inclusive range (start and end coordinates included)
            if orf_length % 3 != 0 and orf not in invalid_orfs[strand]:
                print("Invalid orf: {}; Not multiple of three".format(orf['coords']))
                invalid_orfs[strand].append(orf)
                invalid_msg.append([orf['coords'][0],"Not multiple of three"])

    return invalid_orfs, invalid_msg


def sort_orfs(orf_locations: dict) -> dict:
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
            difference = abs(rev_orf['coords'][0][0] - first_orf['coords'][0][0])  % 3
            if difference == 0:
                sorted_orfs['-0'].append(rev_orf)
            elif difference == 1:
                sorted_orfs['-1'].append(rev_orf)
            elif difference == 2:
                sorted_orfs['-2'].append(rev_orf)

    return sorted_orfs


def set_global_omega_values(orf_locations: dict, omega_values: list, omega_shape: float, omega_classes: int) -> dict:
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

def create_log_file(input_file_name: str) -> str:
    """
    Create a log file with information for the run
    """
    # Select name of the input file
    file_name = input_file_name.split("/")[-1]
    file_name = file_name.split(".")[0]
    LOG_FILENAME = "{}_evol_simulation.log".format(file_name)

    return LOG_FILENAME


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


def count_internal_stop_codons(cds: str) -> int:
    """
    Look for stop codons inside the CDS
    :param cds: str, Coding Aequences (if on negative strand, cds must be the complement)
    :return: the number of stop codons in the coding sequence
    """
    pat = '(TAA|TGA|TAG)'
    reg = re.compile(pat)
    stop_count = 0
    stop_codons = ['TAA', 'TGA', 'TAG']

    # Find if STOP codons in cds by iterating sequence every three nucleotides
    for i in range(3, len(cds)+1, 3):
        codon = cds[(i-3):i]
        if codon in stop_codons:
            stop_count +=1

    return stop_count

def find_ovrfs(orf_list):
    """
    Search for overlapping regions between ORF coordinates
    :param: orf_list: list of list. Start and end positions of reading frames in seq (including spliced and circular CDSs)
    :return: Overlapping regions and ORFs involved 
    """

    # TO DO: How to calculate overlaps on spliced genomes? 
    n = len(orf_list)
    overlaps = {}

    for i in range(n):
        orf1 = orf_list[i]  # Select one orf
        xl1 = min([l for l, r in orf1])  # Find extreme left

        for j in range(i):  # Compare against every other orf
            orf2 = orf_list[j]
            xl2 = min([l for l, r in orf2])

            for l1, r1 in orf1:
                for l2, r2 in orf2:
                    left = max(l1, l2)
                    right = min(r1, r2)
                    overlap = (right - left) +1
                    if overlap > 0:
                        k = str(orf1) + 'and' + str(orf2)
                        overlaps[k] = [orf1, orf2, overlap, xl1, xl2]

    return overlaps

def omegas_in_orf(seq):
        """
        From sequence object, get omega rates in each codon per orf
        :param: Sequence object.
        :return: dict, keyed by orf coordinates. Values are lists with omegas in orf
        """
        codons = seq.get_codons()
        orf_omegas = {}
        for codon in codons:
            orf = str(codon.orf['coords'])
            if orf not in orf_omegas:
                orf_omegas[orf] = [codon.omega]
            else:
                orf_omegas[orf].append(codon.omega)

        return orf_omegas

def main():

    parser = argparse.ArgumentParser(
        description='Simulates and visualizes the evolution of a sequence through a phylogeny'
    )
    args = get_args(parser)
    input = args.seq.lower()

    start_time = datetime.now()
    print("\nStarted at: ", datetime.now())
    
    # Create log file
    file_name = input.split("/")[-1]
    LOG_FILENAME = args.logfile if args.logfile else (f'{file_name.split(".")[0]}_evol_simulation.log')
    logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
    logging.info(f"\nSimulation started at: {start_time}\n")

    # Assign settings
    settings = Settings(args)

    # Read in sequence and resolve ambiguos nucleotides
    s = resolve_ambiguities(str(settings.seq))

    # Get parameters for run
    pi = settings.pi
    global_rate = settings.global_rate
    kappa = settings.kappa
    mu_values = settings.mu_values
    
    pp = pprint.PrettyPrinter(indent=2)
    # Log global parameters
    logging.info(   f"\n\n"
                    f"FILES\n"
                    f"\tSequence: {args.seq}\n"
                    f"\tConfiguration: {args.config}\n"
                    f"\tPhylo Tree: {args.tree}\n"
                    f"\tAlignment: {args.outfile}\n"

                    f"\n"
                    f"PARAMETERS \n"
                    f"\tPi: {pi}\n"
                    f"\tGlobal rate: {global_rate}\n"
                    f"\tKappa: {kappa}\n"
                    f"\tRate classes info: \n{pp.pformat(mu_values)}\n")
    
    orf_locations = settings.orfs
    # Check if the ORFs are valid
    invalid_orfs, invalid_msg = valid_orfs(orf_locations, len(s))

    # Omit the invalid ORFs
    if invalid_orfs['+'] or invalid_orfs['-']:
        logging.warning(f" Omitted ORFS:\n{invalid_msg}\n")
        for strand in invalid_orfs:
            orfs = invalid_orfs[strand]
            for orf in orfs:
                orf_locations[strand].remove(orf)

    # Find internal STOP, remove from orf_locations CDSs with more than 1 STOP codon
    cds = []
    for strand in orf_locations:
        orfs = orf_locations[strand]
        for orf in orfs:
            for start, stop in orf['coords']:
                cds.extend(s[start:stop])  # concatenates spliced ORFs
            
            if strand.startswith('-'):
                cds = cds[::-1]  # negative strand ORF
                cds = Sequence.complement(cds)  # negative strand ORF

            # Count internal STOP codons on CDS
            stop_count = count_internal_stop_codons(cds)
            if stop_count > 1:
                print(f"Omitted orf: {orf['coords']}, has {stop_count} STOP codons")
                orf_locations[strand].remove(orf)  # Remove orf from list of orfs in strand

    # Orf Map
    # Array of one's and cero's used to define position of the orf in a list with as many possitions as orfs in seq
    # It would be latter used as binary code to map a number for the branches on the event tree
    # e.g. [0,1,0] for the second orf in a sequence with three orfs
    orf_coords = [orf['coords'][0] for strand in orf_locations for orf in orf_locations[strand]]
    orf_coords = sorted(orf_coords, key=min)
    orf_maps = np.identity(len(orf_coords), dtype=int)
    orf_map_dict = {tuple(orf_coord):orf_maps[idx,:].copy() for idx, orf_coord in enumerate(orf_coords)}
    for strand in orf_locations:
        for orf in orf_locations[strand]:
            orf['orf_map'] = orf_map_dict[tuple(orf['coords'][0])]

    # Since ORFs are valid, sort the ORFs by reading frame
    orfs = sort_orfs(orf_locations)

    # Final orf list:
    orfs_list = []
    for frame, orf_list in orfs.items():
        for orf in orf_list:
            orfs_list.append(orf['coords'])

    logging.info(f"\n\nVALID ORFs\n"
                 f"{pp.pformat(orfs_list)}\n"
                 f"\nTotal ORFs: {len(orfs_list)}\n"
                 f"\nORFs INFO \n{pp.pformat(orfs)}\n")
    print(f"\nValid orfs: {orfs_list}\nTotal orfs: {len(orfs_list)}\n")
    
    # Check if sequence is valid
    if not valid_sequence(s):
        print("still invalid")
        print("Invalid sequence")
        logging.error("Invalid sequence: {}".format(s))
        sys.exit(0)
   
    # Read in the tree
    phylo_tree = Phylo.read(settings.tree, 'newick', rooted=True)

    # Make Sequence object
    print("Creating root sequence")
    root_sequence = Sequence(s, orfs, kappa, global_rate, pi, mu_values)
    if args.ci:
        omegas = omegas_in_orf(root_sequence)
        for coord, omega_list in omegas.items():
            mean = sum(omega_list)/len(omega_list)
            print(f"coordinates: {coord}, mean omega: {mean}")
        with open(f"{args.outfile}.omegas", 'w+') as write_file:
            json.dump(omegas, write_file, indent=4)

    print(f"Regions info:")
    pp.pprint(root_sequence.regions)
    print("\n")
    logging.info(f"\n\nREGIONS\n{pp.pformat(root_sequence.regions)}\n")

    # Run simulation
    simulation = SimulateOnTree(root_sequence, phylo_tree, args.outfile)
    simulation.get_alignment(args.outfile)

    end_time = datetime.now()
    logging.info(
                    f"\n\tSimulation Ended at: {end_time}\n"
                    f"\tSimulation lasted: {end_time - start_time} seconds\n"
                    f"----------------------------------------------------------\n")

    print(f"Simulation completed.\nAlignment at: {args.outfile}.\nRun Information at:{args.logfile}.\n")
    print(f"Duration: {end_time - start_time} seconds")

if __name__ == '__main__':
    main()
