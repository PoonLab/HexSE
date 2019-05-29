from sequence_info import Sequence
from ovrf_functions import get_frequency_rates
from ovrf_functions import draw_omega_values
from ovrf_functions import get_syn_subs

NUCLEOTIDES = ['A', 'C', 'G', 'T']

class Rates(list):
    """
    List of dictionaries for rates of 3 possible nucleotide substitution at each site in <seq> (3 x sequence_length)
    """

    def __init__(self, seq, mu, bias, pi, omega):
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
    :param nt: nucleotide to which seq[position] will change
    :return: dictionary of updated rates for mutated nucleotide
    """

    rates_before_omega = rates.before_omega[position]
    orfs = rates.seq.orfs
    for i in range(6): #orfs
        print(orfs[i])
        # if self.seq[position].in_orf[i]:
        #     if i < 3:
        #         # positive strand
        #         position in orf
        #

    # #print(rates_before_omega)
    # new_seq = rates.seq[(position - 2):(position + 3)]
    # nt_codons = rates.seq.codon[position]
    # orfs = rates.orfs
    # orfs_for_position = new_seq[2].in_orf
    # #print(orfs, orfs_for_position, "\n", nt_codons, "\n")
    # new_rates = rates[(position - 2):(position + 3)]
    # #print("RATES:", new_rates)
    # #print()






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
