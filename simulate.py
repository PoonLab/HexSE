from Bio import Phylo
import numpy as np
import random
from ovrf_functions import NUCLEOTIDES


class Simulate:
    """
    Simulate evolution within a sequence throughout a phylogeny
    """

    def __init__(self, rates):
        self.rates = rates  # positional list of dictionaries keyed by nt storing rates
        self.total_rate = self.sum_rates()
        self.seq = rates.seq

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

        print(limit)
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

        return(position, to_nt[0])

    def simulate_on_branch(self, evolution_time):
        """
        Simulate molecular evolution on the branch given starting sequence
        """

        # Generate random waiting times to mutate while sum(t)<=branch_length
        times_sum = 0
        while True:
            random_time = random.random()
            times_sum += random_time
            if round(times_sum, 10) > evolution_time:
                break
            # Mutate sequence
            mutation_site = self.get_substitution()[0]
            nucleotide = self.get_substitution()[1]
            self.seq[mutation_site] = nucleotide

        return ''.join(self.seq)

    def traverse_tree(self, tree):
        """
        Post-order traversal of tree from the root.
        Call simulate_on_branch on each branch and feed the resulting sequence to initialize
        the next call.
        @return Phylo tree with Clade objects annotated with sequences.
        """
        # assign root_seq to root Clade
        tree.root.sequence = ''.join(self.seq)
        print(tree)
        # annotate Clades with parents
        for clade in tree.find_clades(order='level'):
            for child in clade:
                child.parent = clade
                print(child, child.parent)

        for node in tree.find_clades(order='level'):
            # TODO: skip the root (sequence already assigned)
            print("--", node)
            if not hasattr(node, 'sequence'):
                # print("second loop", node, node.parent, node.parent.sequence)
                node.sequence = self.simulate_on_branch(node.branch_length)
                print('*', node.sequence)

        # Cleanup to avoid RecursionError when printing the tree
        for clade in tree.find_clades(order='level'):
            for child in clade:
                del child.parent

        return tree

    def get_alignment(self, tree):
        """
        Iterates over tips (terminal nodes) of tree and returns sequence
        attributes as list.
        :param tree:
        :return:
        """
