import scipy
from scipy.stats import gamma
import numpy as np

# Calculate evolution rates across the sequence
def get_frequency_rates(seq):
    seq = list(seq)
    frequencies = {
        'A': 0,
        'C': 0,
        'T': 0,
        'G': 0,
    }

    for nucleotide in frequencies:
        frequencies[nucleotide] = round((float(seq.count(nucleotide))/(len(seq))),2)
    return frequencies


def evol_rates(seq, mu, bias, pi, orfs):
    """
    Generate rate vector from sequence given parameters.
    @param seq: the nucleotide sequence of length <L>.  Assumed to start and end in reading frame +0.
    @param mu: the global rate (substitutions/site/unit time).
    @param pi: <optional> a vector of stationary nucleotide frequencies.  If not specified, defaults
               to the empirical frequencies in <seq>.
    @param bias: <optional> a List of 6 rates [AC, AG, AT, CG, CT, GT] assuming time-reversibility.
                    If not specified, defaults to 1's.
    @param omega: <optional> a List of dN/dS (omega) ratios along sequence length as {dict}.
                  {dict} always contains key +0 (parent reading frame).  May contain omega for alternate
                  reading frame as (+1, +2, -0, -1 or -2).  Codon position is determined by the nt's
                  location relative to start of <seq>.
    @param orf: tuple indicated by user containing first and last nucleotide of every reading frame in seq (ex. ((4,24),(6,17))
    @return seq_rates: a List of tuples for rates of 3 possible nucleotide substitution at each site
                   in <seq> (3xL).
    """

    L = len(seq)
    seq_rates = [(mu, mu, mu) for position in range(L)]                      # initialize the list with baseline rates
    frequency_rates = get_frequency_rates(seq)
    # iterate over every nt in seq
    for position in range(L):
        # 1. what is the current nucleotide X?
        nucleotide = seq[position]

        # 2. apply stationary frequencies (pi) to tuple given X
        if pi is None:                                                       # If pi is no specified
            pi = frequency_rates[nucleotide]                                 # use empirical frequencies

        seq_rates[position] = tuple([pi*j for j in seq_rates[position]])     # Apply pi to every value in seq_rates according to the position

        # 3. apply biases to tuple
        biased_rates = []

        for rate in range(3):                                                      # Iterate over every item in the tuple
            modified_rates = seq_rates[position][rate] * bias[nucleotide][rate]    # Create a new list with values of the tuple modified by bias
            biased_rates.append(modified_rates)

        seq_rates[position] = tuple(biased_rates)                        # Update seq_rates

        # 4. Apply omega

        #Check. If for each codon if substitution is nonsynonymous. If true, apply omega (over the codon or the nt).
    return seq_rates


def get_omega(orf):
    """
    Draw omega values for every reading frame in seq from a gamma distribution
    @param orf: tuple indicated by user containing first and last nucleotide of every reading frame in seq (ex. ((4,24),(6,17))
    @return omega: dictionary with keys as beginning and end of the RF in seq and the dN/dS rates for each codon.
     """
    omega = {}
    a = 1                       # Shape parameter
    for i in orf:
        number_of_codons = (i[1] - i[0])//3
        omega[i] = gamma.rvs(a, size = number_of_codons)
    return omega


def get_reading_frames(seq):
    """
    Creates a list with tuples containing the first and last position of the reading frames in seq
    according to start and stop codons
    """
    start = 'ATG'
    stop = ['TAG', 'TAA']
    reading_frames = []
    for frame in range(3):
        for codon, position in codon_iterator(seq[frame:]): # Iterate over every codon in the RF
            if codon == start:                              # Find a start codon
                for codon_in_rf, position_in_rf in codon_iterator(seq[(position+frame):]):
                    if codon_in_rf in stop:                 # Find a stop codon
                        # Get the positions in the sequence for the first and last nt of the RF
                        orf = (position+frame, position_in_rf+(position+frame)+3)
                        #Use +3 to include the full stop codon
                        reading_frames.append(orf)
                        break
                break
    return reading_frames


def codon_iterator(seq_in_orf):
    """
    Generator to move every tree nucleotides (codon)
    @param seq_in_orf: nucleotide sequence in a given orf
    Yield: codon and position of the first nucleotide of codon in seq (ex., ('GAA', 5))
    """
    i = 0
    while i < len(seq_in_orf):
        yield list[i:i + 3], i
        i += 3

def get_syn_codons(my_codon):
    """
    Create a list with synonymous codons of a given codon
    @param codon: Three nucleotides in an ORF
    """
    my_aa = codon_dict[my_codon]
    syn_codons = []
    for codon, aa in codon_dict.items():
        if codon_dict[codon] == my_aa:
            syn_codons.append(codon)
    return syn_codons

# Problem: it is only accounting for plus reading frames
def get_syn_subs(seq, orfs, codon_dict):
    """
    Get synonynous and non-synonymous substitutions for nucleotide in seq, given the open reading frames
    @param seq: nucleotide sequence
    @param orfs: tuple of open reading frames present in the sequence (ex: ([0,11],[1,6]))
    @param codon_dict: dictionary with translated codons
    @return: list of dictionaries with possible mutations. If mutation is syn, then value for that mutation is zero.
             If the mutation is nonsyn, value is one.
    """
    L = len(seq)
    nts = ['A', 'C', 'G', 'T']
    seq_subs = []                  # Possible substitution for the entire sequence
    for position in range(L):
        nt_subs = {}               # Possible substitutions for the current nucleotide
        for i in range(len(nts)):
            is_syn = []            # Checking if the substitution is syn or nonsyn
            to_nt = nts[i]
            for j in range(len(orfs)):
                orf = orfs[j]
                if position in range(orf[0],orf[1]):  # is nt in orf?
                    nt_info = get_codon(seq, position, orf)
                    my_codon = nt_info[0]
                    position_in_codon = nt_info[1]
                    mutated_codon = list(my_codon)
                    mutated_codon[position_in_codon] = to_nt # substitution step
                    if codon_dict[''.join(mutated_codon)] == codon_dict[my_codon]:
                        is_syn.append(0)
                    else:
                        is_syn.append(1)
                else:
                    is_syn.append(0)

            nt_subs[to_nt] = is_syn

        seq_subs.append(nt_subs)

    return seq_subs


def get_codon(seq, position , orf):
    """
    Get codon, and position in that codon of a nucleotide given an orf
    @param seq: parental sequence
    @param position: position of the nucleotide in <seq>
    @param orf: tuple indicated by user containing first and last nucleotide of an open reading frame
    @return position_in_codon: of the current nucleotide in the triplet
    @return codon: to which the nucleotide belongs to
    """

    my_orf = seq[orf[0]:orf[1]+1]
    position_in_orf = position - orf[0]
    # if position_in_orf < 0, then raise argument error
    if position_in_orf % 3 == 0:
        position_in_codon = 0
        codon = my_orf[position_in_orf:position_in_orf+3]
    elif position_in_orf %3 == 1:
        position_in_codon = 1
        codon = my_orf[position_in_orf-1:position_in_orf+2]
    else:
        position_in_codon = 2
        codon = my_orf[position_in_orf-2:position_in_orf+1]

    return codon, position_in_codon


