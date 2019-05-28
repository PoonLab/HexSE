import scipy
from scipy.stats import gamma
import numpy as np
from sequence_info import Sequence

COMPLEMENT_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'W', 'R': 'Y', 'K': 'M', 'Y': 'R',
                   'S': 'S', 'M': 'K', 'B': 'V', 'D': 'H',
                   'H': 'D', 'V': 'B', '*': '*', 'N': 'N',
                   '-': '-'}

NUCLEOTIDES = ['A', 'C', 'G', 'T']

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


class Rates(list):
    """
    List of dictionaries for rates of 3 possible nucleotide substitution at each site in <seq> (3 x sequence_length)
    """

    def __init__(self, seq, mu, sorted_orfs, bias, pi, omega):
        """
        Generate rate vector from sequence given parameters.
        :param seq: the nucleotide.
        :param mu: the global rate (substitutions/site/unit time).
        :param orfs: list of tuples indicated by user containing first and last nucleotide
                        of every reading frame in seq (ex. [(4,24),(6,17]) if there is not orf in some possible reading frame, it will have a "1"
        :param bias: <optional> a List of 6 rates [AC, AG, AT, CG, CT, GT] assuming time-reversibility.
                        If not specified, defaults to 1's.
        :param pi: <optional> a vector of stationary nucleotide frequencies.  If not specified, defaults
                   to the empirical frequencies in <seq>.
        :param omega: <optional> a List of dN/dS (omega) ratios along sequence length as {dict}.
                      {dict} always contains key +0 (parent reading frame).  May contain omega for alternate
                      reading frame as (+1, +2, -0, -1 or -2). Codon position is determined by the nt's
                      location relative to start of <seq>.
        """

        super(Rates, self).__init__()
        self.seq = Sequence(seq, sorted_orfs)
        self.mu = mu
        self.bias = bias
        self.pi = pi
        self.orfs = sorted_orfs
        self.omega = omega

        #print(self.seq, self.mu, self.bias, self.pi, self.orfs, "\n")

        # Defining parameters:
        if self.pi is None:
            self.pi = get_frequency_rates(self.seq)

        if self.omega is None: # user does not specify omega values
            self.omega = draw_omega_values(self.orfs) # omega values for each orf drawn from gamma distribution

        syn_nonsyn = get_syn_subs(self.seq, self.orfs)

        for position in range(len(self.seq)):
            my_nt = self.seq[position]
            sub_rate = {} # Substitution rates

            # Apply mu, bias an pi
            for to_nt in NUCLEOTIDES:
                if to_nt == my_nt:
                    sub_rate[to_nt] = None
                else:
                    sub_rate[to_nt] = self.mu * self.bias[my_nt][to_nt] * self.pi[my_nt]
            self.append(sub_rate)

            # Apply omega
            for i in range(len(self.orfs)):
                orf = self.orfs[i]
                if type(orf) == tuple:
                    omega_values = self.omega[orf]

                    if position in range(orf[0], orf[1] + 1):  # if nt in orf (positive strand)
                        position_in_orf = position - orf[0]
                        codon_in_orf = position_in_orf // 3
                        for to_nt in NUCLEOTIDES:
                            if self[position][to_nt] is not None:  # There is a reading frame
                                if syn_nonsyn[position][to_nt][i] == 1:  # mutation is Non-synonymous
                                    self[position][to_nt] *= omega_values[codon_in_orf]  # Apply omega

                    elif position in range(orf[0] + 1, orf[1], -1): # if nt in orf (negative strand)
                        position_in_orf = orf[0] - position
                        codon_in_orf = position_in_orf // 3
                        for to_nt in NUCLEOTIDES:
                            if self[position][to_nt] is not None:  # There is a reading frame
                                if syn_nonsyn[position][to_nt][i] == 1:  # mutation is Non-synonymous
                                    self[position][to_nt] *= omega_values[codon_in_orf]  # Apply omega


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
        'G': 0,
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
        if type(i) == tuple: # from sorted_orfs, if there is no orf in some position, skip it
            if i[1] > i[0]:
                number_of_codons = (((i[1] + 1) - i[0]) // 3)
            else:
                number_of_codons = (((i[0] + 1) - i[1]) // 3)

            omega_values[i] = gamma.rvs(a, size = number_of_codons)

    return omega_values


def get_reading_frames(seq):
    """
    Creates a list with tuples containing the first and last position
    of the reading frames in seq according to START and STOP codons
    :param seq: the input DNA sequence
    :return: a list of tuples containing the index of the first nucleotide of
                the START codon and the index of the last nucleotide of the STOP codon
    """
    start = 'ATG'
    stop = ['TAG', 'TAA', 'TGA']
    reading_frames = []
    # Make an if in case that no stop or start are found
    for frame in range(3):
        for codon, position in codon_iterator(seq[frame:]):  # Iterate over every codon in the RF
            if codon.upper() == start:  # Find a start codon
                for codon_in_rf, position_in_rf in codon_iterator(seq[(position + frame):]):
                    if codon_in_rf.upper() in stop:  # Find a stop codon
                        # Get the positions in the sequence for the first and last nt of the RF
                        orf = (position + frame, position_in_rf + (position + frame) + 2)
                        # Use +3 to include the full stop codon
                        reading_frames.append(orf)
                        break
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
    my_aa = CODON_DICT[my_codon.upper()]
    syn_codons = []
    for codon, aa in CODON_DICT.items():
        if CODON_DICT[codon] == my_aa:
            syn_codons.append(codon)

    return syn_codons


def get_syn_subs(seq, orfs):
    """
    Get synonynous and non-synonymous substitutions for nucleotide in seq, given the open reading frames
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
                if type(orf) == tuple: # does it exist?
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
                else: # there is no reading frame
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
    if orf[1] > orf[0]: # positive strand
        my_orf = ''.join(seq[orf[0]:orf[1] + 1])
        position_in_orf = position - orf[0]
    else: # negative strand
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


def update_rates(rates, position, nt):
    """
    :param rates: instance of Rates for <seq>.
    :param position: position where the mutation occurs
    :param nt: nucleotide to which seq[position] will change
    :return: dictionary of updated rates for mutated nucleotide
    """

    new_seq = rates.seq[(position - 2):(position + 3)]
    orfs = rates.orfs
    orfs_for_position = new_seq[2].in_orf
    print(orfs, orfs_for_position)
    new_rates = rates[(position - 2):(position + 3)]
    print(new_rates)






    # for orf in orfs:
    #     temp = get_codon(seq, position, orf) # get position and codon in which seq[position] is involved
    #     codon = list((temp)[0])
    #     pos_in_codon = temp[1]
    #     if orf[1] > orf[0]: # Positive strand
    #         to_nt = nt
    #     else:
    #         to_nt = complement_dict[nt]
    #
    # nt_rates = Rates(mutated_codon, rates.mu, rates.bias, rates.pi, local_orfs)[pos_in_codon] # get new rates
    # return(nt_rates)
    #
    #
