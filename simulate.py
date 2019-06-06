from Bio import Phylo
import numpy as np
import random
from ovrf_functions import NUCLEOTIDES
from ovrf_functions import COMPLEMENT_DICT
from ovrf_functions import CODON_DICT


class Simulate:
    """
    Simulate evolution within a sequence throughout a phylogeny
    """

    def __init__(self, rates, tree):
        self.rates = rates  # positional list of dictionaries keyed by nt storing rates
        self.total_rate = self.sum_rates()
        self.seq = rates.seq
        self.tree = tree

    def sum_rates(self):
        """
        Calculate the total rate by iterating over dictionaries in <rates>
        """
        res = 0
        for pdict in self.rates:
            for nt, rate in pdict.items():
                if rate != None:
                    res += rate
        return res

    def get_substitution(self):
        """
        Randomly select a position, and a nucleotide where a mutation in <seq> occurs, given substitution rates
        @return: position and nucleotide to which <seq> will change
        """
        # Draw a random limit where to reach te mutation given rates
        limit = random.uniform(0, self.sum_rates())
        position = -1 # Should include position zero and can only go to len(seq)-1
        total = 0
        # draw position
        while total < limit:
            pdict = self.rates[position]
            for to_nt, rate in pdict.items():
                if rate is not None:
                    total += rate
            position += 1


        # draw nucleotide
        local_rates = self.rates[position]
        nucleotides = []
        probabilities = []
        for nt, rate in local_rates.items():
            if rate is not None:
                nucleotides.append(nt)
                probabilities.append(rate / self.total_rate)

        to_nt = np.random.choice(NUCLEOTIDES, 1, probabilities)

        # update rates if position is in at least in one reading frame
        total = []
        for orf in self.rates.orfs:
            if type(orf) == tuple:
                total.append(orf[1])
                total.append(orf[0])

        if position in range(min(total), max(total)):
            self.rates[position] = self.get_new_rates(position, to_nt[0])

        return(position, to_nt[0])

    def simulate_on_branch(self, evolution_time):
        """
        Simulate molecular evolution on the branch given starting sequence
        """

        # Generate random waiting times to mutate while sum(t)<=branch_length
        times_sum = 0
        while True:
            random_time = random.uniform(0, evolution_time)
            #random_time = np.random.exponential(scale=1, size = None)
            times_sum += random_time
            if round(times_sum, 10) > evolution_time:
                break
            # Mutate sequence
            mutation_site = self.get_substitution()[0]
            nucleotide = self.get_substitution()[1]
            self.seq[mutation_site] = nucleotide

        return ''.join(self.seq)

    def traverse_tree(self):
        """
        Post-order traversal of tree from the root.
        Call simulate_on_branch on each branch and feed the resulting sequence to initialize
        the next call.
        @return Phylo tree with Clade objects annotated with sequences.
        """
        # assign root_seq to root Clade
        self.tree.root.sequence = ''.join(self.seq)

        # annotate Clades with parents
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                child.parent = clade
                #print(child, child.parent)

        for node in self.tree.find_clades(order='level'):
            # TODO: skip the root (sequence already assigned)
            # print("--", node)
            if not hasattr(node, 'sequence'):
                # print("second loop", node, node.parent, node.parent.sequence)
                node.sequence = self.simulate_on_branch(node.branch_length)
                #print('*', node.sequence)

        # Cleanup to avoid RecursionError when printing the tree
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                del child.parent

        self.final_tree = self.tree
        return self.final_tree

    def get_alignment(self):
        """
        Iterates over tips (terminal nodes) of tree and returns sequence
        attributes as list.
        :param tree:
        :return:
        """
        for node in self.final_tree.find_clades(order='level'):
            print(node.sequence)


    def get_new_rates(self, position, nt):
        """
        Recalculate substitution rates for the mutated nucleotide
        :param rates: instance of Rates for <seq>.
        :param position: position where the mutation occurs
        :param nt: nucleotide to which seq[position] will change (mutation)
        :return: dictionary of updated rates for mutated nucleotide
        """

        print("Running get_reading_frames")
        print(position, nt)
        # print(self.rates[position])
        rates_before_omega = self.rates.before_omega[position]
        all_codons_info = self.rates.seq.codon[position]  # Every codon in which nucleotide is involved given ORFs
        # print("this are the codons")
        # print(all_codons_info)
        # print("this are the orfs")
        # print(self.rates.orfs)
        new_rates = {}

        # Change codons in which nt is involved for new mutated nucleotide
        mutated_codons_info = []
        for i in range(6):

            if all_codons_info[i] != 0:
                # if there is a codon for the orf, extract codon information
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
                    # print(rates.mu, rates.bias[nt][to_nt], rates.pi[nt])
                    rates_before_omega[to_nt] = self.rates.mu * self.rates.bias[nt][to_nt] * self.rates.pi[nt]

        rates_before_omega[nt] = None

        # Apply omega
        for i in range(6):
            specific_orf = self.rates.orfs[i]

            if type(specific_orf) == tuple:
                # means that this is an ORF, otherwise no shift in the seq
                omega_values = self.rates.omega[specific_orf]

                if self.rates.seq[position].in_orf[i]: # if nt in orf
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
                                # Is a non-synonymous mutation
                                new_rates[to_nt] = rates_before_omega[to_nt] * omega_values[codon_in_orf]

                            else:
                                new_rates[to_nt] = rates_before_omega[to_nt]

                        else:
                            new_rates[to_nt] = None

        # Update mutated codons
        self.seq.codon[position] = mutated_codons_info

        return new_rates


