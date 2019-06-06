from sequence_info import Sequence
from ovrf_functions import NUCLEOTIDES
from ovrf_functions import COMPLEMENT_DICT
from ovrf_functions import CODON_DICT
from scipy.stats import gamma


class Rates(list):
    """
    List of dictionaries for rates of 3 possible nucleotide substitution at each site in <seq> (3 x sequence_length)
    """

    def __init__(self, seq, mu, bias=(1, 1, 1, 1, 1, 1), pi=None, omega=None):
        """
        Generate rate vector from sequence given parameters.
        :param seq: the nucleotide.
        :param mu: the global rate (substitutions/site/unit time).
        :param orfs: list of tuples indicated by user containing first and last nucleotide
                        of every reading frame in seq (ex. [(4,24),(6,17])
                        if there is not orf in some possible reading frame, it will have a "1"
        :param bias: <optional> a List of 6 rates [AC, AG, AT, CG, CT, GT] assuming time-reversibility.
                        If not specified, defaults to 1's.
        :param pi: <optional> a vector of stationary nucleotide frequencies.  If not specified, defaults
                   to the empirical frequencies in <seq>.
        :param omega: <optional> a List of dN/dS (omega) ratios along sequence length as {dict}.
                      {dict} always contains key +0 (parent reading frame).  May contain omega for alternate
                      reading frame as (+1, +2, -0,sim_ -1 or -2). Codon position is determined by the nt's
                      location relative to start of <seq>.
        """

        super(Rates, self).__init__()
        self.seq = Sequence(seq, seq.orfs)
        self.mu = mu
        self.bias = bias
        self.pi = pi
        self.orfs = self.seq.orfs
        self.omega = omega

        #print(self.seq, self.mu, self.bias, self.pi, self.orfs, "\n")

        # Defining parameters:
        if self.pi is None:
            self.pi = self.get_frequency_rates(self.seq)

        if self.omega is None: # user does not specify omega values
            self.omega = self.draw_omega_values(self.orfs) # omega values for each orf drawn from gamma distribution

        syn_nonsyn = self.get_syn_subs(self.seq, self.orfs)

        self.before_omega = []
        for position in range(len(self.seq)):
            my_nt = self.seq[position]
            sub_rate = {} # Substitution rates

            # Apply mu, bias an pi
            for to_nt in NUCLEOTIDES: #create dict
                if to_nt == my_nt:
                    sub_rate[to_nt] = None
                else:
                    sub_rate[to_nt] = self.mu * self.bias[my_nt][to_nt] * self.pi[my_nt]

            self.before_omega.append(sub_rate.copy())
            self.append(sub_rate)

            # Apply omega
            for i in range(len(self.orfs)):
                specific_orf = self.orfs[i]

                if type(specific_orf) == tuple:
                    # means that this is an ORF, otherwise no shift in the seq
                    omega_values = self.omega[specific_orf]

                    if self.seq[position].in_orf[i]:
                        if i < 3:
                            # positive strand
                            position_in_orf = position - specific_orf[0]
                        else:
                            # negative strand
                            position_in_orf = specific_orf[0] - position

                        codon_in_orf = position_in_orf // 3
                        for to_nt in NUCLEOTIDES:
                            if self[position][to_nt] is not None:
                                # access the rate, and modify if non-synonymous
                                if syn_nonsyn[position][to_nt][i] == 1:
                                    self[position][to_nt] *= omega_values[codon_in_orf]

    def draw_omega_values(self, orfs):
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

    def get_syn_subs(self, seq, orfs):
        """
        Get synonymous and non-synonymous substitutions for nucleotide in seq, given the open reading frames
        :param seq: nucleotide sequence
        :param orfs: tuple of open reading frames present in the sequence (ex: [(0,11),(1, 6)])
        :return: list of dictionaries with possible mutations. If mutation is syn, then value for that mutation is zero.
                 If the mutation is nonsyn, value is one.
        """

        seq_subs = []  # Possible substitution for the entire sequence

        for pos, nucl in enumerate(seq):
            nt_subs = {}
            # Possible substitutions for the current nucleotide
            for nt in NUCLEOTIDES:
                is_syn = []  # Checking if the substitution is syn or nonsyn
                for i in range(6):
                    orf = orfs[i]
                    if type(orf) == tuple: # if orf exists
                        if orf[0] <= pos <= orf[1] or orf[0] >= pos >= orf[1]: # if position in orf
                            nt_info = seq.codon[pos][i]
                            my_codon = nt_info[0]
                            position_in_codon = nt_info[1]
                            mutated_codon = list(my_codon)
                            mutated_codon[position_in_codon] = nt  # substitution step

                            is_syn.append(0) if CODON_DICT[''.join(mutated_codon)] == CODON_DICT[my_codon] \
                                else is_syn.append(1)
                        else:
                            is_syn.append(0)
                    else:
                        is_syn.append(0)

                nt_subs[nt] = is_syn
            seq_subs.append(nt_subs)

        return seq_subs

    def get_frequency_rates(self, seq):
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

