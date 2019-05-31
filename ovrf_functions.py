# Script with the necessary functions to calculate evolution rates and extract sequence information
# Prepare input for the simulation script

import re

import numpy as np
import scipy
from scipy.stats import gamma

COMPLEMENT_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'W', 'R': 'Y', 'K': 'M', 'Y': 'R',
                   'S': 'S', 'M': 'K', 'B': 'V', 'D': 'H',
                   'H': 'D', 'V': 'B', '*': '*', 'N': 'N',
                   '-': '-'}

CODON_DICT = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
              'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
              'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
              'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
              'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
              'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
              'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
              '---': '-', 'XXX': '?'}

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def get_frequency_rates(seq):
    """
    Frequency of nucleotides in the DNA sequence
    :param seq: the DNA sequence
    :return: a dictionary frequencies where the key is the nucleotide and the value is the frequency
    """
    seq = list(seq)
    frequencies = {
        'A': 0,
        'C': 0,
        'T': 0,
        'G': 0
    }

    for nucleotide in frequencies:
        frequencies[nucleotide] = round((float(seq.count(nucleotide)) / (len(seq))), 2)

    return frequencies


def draw_omega_values(orfs):
    """
    Draw omega values for every reading frame in seq from a gamma distribution
    :param orfs: List of tuples indicated by user containing first and last nucleotide of
                every reading frame in <seq> (ex. [(4,24), (3,15)])
    :return omega: dictionary with keys as beginning and end of the RF in seq and the dN/dS rates for each codon.
     """
    omega_values = {}
    a = 1  # Shape parameter
    for i in orfs:
        if type(i) == tuple:  # from sorted_orfs, if there is no orf in some position, skip it
            if i[1] > i[0]:
                number_of_codons = (((i[1] + 1) - i[0]) // 3)
            else:
                number_of_codons = (((i[0] + 1) - i[1]) // 3)

            omega_values[i] = gamma.rvs(a, size=number_of_codons)

    return omega_values


def get_reading_frames(seq):
    """
    Creates a list with tuples containing the first and last position
    of the reading frames in seq according to START and STOP codons
    :param seq: the input DNA sequence
    :return: a list of tuples containing the index of the first nucleotide of
                the START codon and the index of the last nucleotide of the STOP codon
    """
    start_codon = re.compile('ATG', flags=re.IGNORECASE)
    stop = re.compile('(TAG)|(TAA)|(TGA)', flags=re.IGNORECASE)
    reading_frames = []

    # Record positions of all potential START codons
    start_positions = [match.start() for match in start_codon.finditer(seq)]

    for position in start_positions:
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
                # Find a stop codon and ensure ORF length is sufficient
                if match.start() % 3 == frame and orf_length >= 8:
                    # Get the positions in the sequence for the first and last nt of the RF
                    orf = (position, match.end() - 1)
                    reading_frames.append(orf)
                    break

    return reading_frames


def codon_iterator(seq_in_orf):
    """
    Generator to move every tree nucleotides (codon)
    :param seq_in_orf: nucleotide sequence in a given orf
    Yield: codon and position of the first nucleotide of codon in seq (ex., ('GAA', 5))
    """
    i = 0
    while i < len(seq_in_orf):
        yield seq_in_orf[i:i + 3], i
        i += 3


def get_syn_codons(my_codon):
    """
    Create a list with synonymous codons of a given codon
    :param my_codon: Three nucleotides in an ORF
    :return: a list of synonomous codons
    """
    try:
        my_aa = CODON_DICT[my_codon.upper()]
    except KeyError as e:
        raise KeyError("Invalid codon {}".format(my_codon))

    syn_codons = []
    for codon, aa in CODON_DICT.items():
        if CODON_DICT[codon] == my_aa:
            syn_codons.append(codon)

    return syn_codons


def get_syn_subs(seq, orfs):
    """
    Get synonymous and non-synonymous substitutions for nucleotide in seq, given the open reading frames
    :param seq: nucleotide sequence
    :param orfs: tuple of open reading frames present in the sequence (ex: ([0,11],[1,6]))
    :return: list of dictionaries with possible mutations. If mutation is syn, then value for that mutation is zero.
             If the mutation is nonsyn, value is one.
    """
    sequence_length = len(seq)
    nts = ['A', 'C', 'G', 'T']
    seq_subs = []  # Possible substitution for the entire sequence
    for position in range(sequence_length):
        nt_subs = {}  # Possible substitutions for the current nucleotide
        for i in range(len(nts)):
            is_syn = []  # Checking if the substitution is syn or nonsyn
            to_nt = nts[i]
            for j in range(len(orfs)):
                orf = orfs[j]
                if type(orf) == tuple:  # does it exist?
                    if position in range(orf[0], orf[1]) or position in range(orf[0], orf[1], -1):  # is nt in orf?
                        nt_info = get_codon(seq, position, orf)
                        my_codon = nt_info[0]
                        position_in_codon = nt_info[1]
                        mutated_codon = list(my_codon)
                        mutated_codon[position_in_codon] = to_nt  # substitution step
                        if CODON_DICT[''.join(mutated_codon)] == CODON_DICT[my_codon]:
                            is_syn.append(0)
                        else:
                            is_syn.append(1)
                    else:
                        is_syn.append(0)
                else:  # there is no reading frame
                    is_syn.append(0)

            nt_subs[to_nt] = is_syn

        seq_subs.append(nt_subs)

    return seq_subs


def get_codon(seq, position, orf):
    """
    Get codon sequence, and position of my_nt in the codon
    :param seq: parental sequence as list of nucleotides
    :param position: position of the nucleotide in <seq>
    :param orf: tuple indicating first and last nucleotide of an open reading frame
    :return codon: tuple with nucleotide triplet and position of the nucleotide in the codon
    """
    if orf[1] > orf[0]:  # positive strand
        my_orf = ''.join(seq[orf[0]:orf[1] + 1])
        position_in_orf = position - orf[0]
    else:  # negative strand
        rseq = reverse_and_complement(seq)
        my_orf = rseq[orf[1]:orf[0] + 1]
        position_in_orf = orf[0] - position

    # if position_in_orf < 0, then raise argument error
    if position_in_orf % 3 == 0:
        position_in_codon = 0
        codon = my_orf[position_in_orf:position_in_orf + 3]
    elif position_in_orf % 3 == 1:
        position_in_codon = 1
        codon = my_orf[position_in_orf - 1:position_in_orf + 2]
    else:
        position_in_codon = 2
        codon = my_orf[position_in_orf - 2:position_in_orf + 1]

    # print (codon, position_in_codon, position_in_orf, "\n")
    return codon, position_in_codon


def reverse_and_complement(seq):
    """
    Generates the reverse complement of a DNA sequence
    :param seq: the DNA sequence
    :return: the reverse complement of the sequence
    """
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        try:
            rcseq += COMPLEMENT_DICT[i]
        except KeyError as e:
            raise KeyError("Invalid character '{}' in sequence".format(i))

    return rcseq


def sort_orfs(orfs):
    """
    Store orfs in position according to plus zero orf (first of the list).
    They will be classified as (+0, +1, +2, -0, -1, -2)
    :param orfs: list of orfs as tuples for <seq> (ex. [(5,16),(11,0)])
    :return: list of ORFs classified according to their shift regarding to the plus zero one (+0, +1, +2, -0, -1, -2)
    """
    plus_cero_orf = orfs[0]
    orf_position = [1] * 6
    orf_position[0] = plus_cero_orf

    for orf in orfs[1:]:
        if type(orf) == tuple:
            if orf[0] < orf[1]:  # positive strand
                difference = (orf[0] - plus_cero_orf[0]) % 3
                if difference == 1:  # plus one
                    orf_position[1] = orf
                elif difference == 2:  # plus two
                    orf_position[2] = orf
            elif orf[0] > orf[1]:  # negative strand
                difference = (plus_cero_orf[1] - orf[0]) % 3
                if difference == 0:
                    orf_position[3] = orf  # minus zero
                elif difference == 1:
                    orf_position[4] = orf  # minus one
                elif difference == 2:
                    orf_position[5] = orf  # minus two

    return orf_position
