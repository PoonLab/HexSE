from sequence_info import Sequence
from ovrf_functions import get_frequency_rates
from ovrf_functions import draw_omega_values
from ovrf_functions import get_syn_subs
from ovrf_functions import NUCLEOTIDES
from ovrf_functions import COMPLEMENT_DICT
from ovrf_functions import CODON_DICT



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
                        of every reading frame in seq (ex. [(4,24),(6,17]) if there is not orf in some possible reading frame, it will have a "1"
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
            self.pi = get_frequency_rates(self.seq)

        if self.omega is None: # user does not specify omega values
            self.omega = draw_omega_values(self.orfs) # omega values for each orf drawn from gamma distribution

        syn_nonsyn = get_syn_subs(self.seq, self.orfs)

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



def update_rates(rates, position, nt):
    """
    :param rates: instance of Rates for <seq>.
    :param position: position where the mutation occurs
    :param nt: nucleotide to which seq[position] will change (mutation)
    :return: dictionary of updated rates for mutated nucleotide
    """

    rates_before_omega = rates.before_omega[position]
    new_rates = {}
    all_codons_info = rates.seq.codon[position]
    omega = rates.omega
    orfs = rates.seq.orfs

    # Change codons in which nt is involved for new mutated nucleotide
    mutated_codons_info = []
    for i in range(6):

        if all_codons_info[i] != 0:
            # extract codon information
            codon = list(all_codons_info[i][0])
            position_in_codon = all_codons_info[i][1]

            if i < 3:
                # positive strand
                codon[position_in_codon] = nt
            else:
                # negative strand
                mutated_nt = COMPLEMENT_DICT[nt]
                codon[position_in_codon] = mutated_nt

            mutated_codon = ''.join(codon)
            new_codon_info = (mutated_codon, position_in_codon)

        else:
            new_codon_info = 0

        mutated_codons_info.append(new_codon_info)

        for to_nt in NUCLEOTIDES:
            if rates_before_omega[to_nt] == None:
                #print(rates.mu, rates.bias[nt][to_nt], rates.pi[nt])
                rates_before_omega[to_nt] = rates.mu * rates.bias[nt][to_nt] * rates.pi[nt]

    rates_before_omega[nt] = None

    # Apply omega

    for i in range(6):
        specific_orf = rates.orfs[i]

        if type(specific_orf) == tuple:
            # means that this is an ORF, otherwise no shift in the seq
            omega_values = rates.omega[specific_orf]

            if rates.seq[position].in_orf[i]:
                if i < 3:
                    # positive strand
                    position_in_orf = position - specific_orf[0]
                else:
                    # negative strand
                    position_in_orf = specific_orf[0] - position

                codon_in_orf = position_in_orf // 3

                for to_nt in NUCLEOTIDES:
                    if rates_before_omega[to_nt] is not None:
                        # access the rate to nt different to itself

                        codon_for_rates = list(mutated_codons_info[i][0])
                        position_in_codon_rates = mutated_codons_info[i][1]

                        codon_for_rates[position_in_codon_rates] = to_nt
                        string_codon = ''.join(codon_for_rates)

                        if CODON_DICT[mutated_codons_info[i][0]] == CODON_DICT[string_codon]:
                            # Is a non-synonymoys mutation
                            new_rates[to_nt] = rates_before_omega[to_nt] * omega_values[codon_in_orf]

                        else:
                            new_rates[to_nt] = rates_before_omega[to_nt]

                    else:
                        new_rates[to_nt] = None