import unittest
import random
import os

from Bio import Phylo
from refact.new_sequence_info import Sequence
from refact.simulation_refactore import SimulateOnBranch
from refact.simulation_refactore import SimulateOnTree

TEST_TREE = os.path.abspath('../refact/test_tree.txt')


# ==========================================
# Tests for SimulateOnBranch
# ==========================================
class TestSimulateOnBranch(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        random.seed(9001)  # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas)
        branch_length = 10
        self.sim_on_branch1 = SimulateOnBranch(sequence1, branch_length)

        s2 = 'TTTTTTCTTTTTTT'
        pi2 = Sequence.get_frequency_rates(s2)
        sequence2 = Sequence(s2, {'+0': [(0, 11)]}, kappa, mu, pi2, omegas)
        self.sim_on_branch2 = SimulateOnBranch(sequence2, branch_length)

        s3 = 'AATTCATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCTAGGCCTACCC'
        sorted_orfs = {'+0': [(5, 49)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': [(29, 3)]}
        pi3 = Sequence.get_frequency_rates(s3)
        sequence3 = Sequence(s3, sorted_orfs, kappa, mu, pi3, omegas)
        self.sim_on_branch3 = SimulateOnBranch(sequence3, branch_length)

        s4 = 'ATGACGTGGTGA'
        sorted_orfs = {'+0': [(0, 11)], '+1': [], '+2': [], '-0': [], '-1': [], '-2': []}
        pi4 = Sequence.get_frequency_rates(s4)
        sequence4 = Sequence(s4, sorted_orfs, kappa, mu, pi4, omegas)
        self.sim_on_branch4 = SimulateOnBranch(sequence4, branch_length)

    def testGetSubstitution(self):
        random.seed(9001)       # Set seed for pseudo-random number generator

        expected = (self.sim_on_branch1.sequence.nt_sequence.nucleotide_at_pos(7), 'A')
        result = self.sim_on_branch1.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch2.sequence.nt_sequence.nucleotide_at_pos(0), 'C')
        result = self.sim_on_branch2.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch3.sequence.nt_sequence.nucleotide_at_pos(6), 'C')
        result = self.sim_on_branch3.get_substitution()
        self.assertEqual(expected, result)

        expected = (self.sim_on_branch4.sequence.nt_sequence.nucleotide_at_pos(8), 'A')
        result = self.sim_on_branch4.get_substitution()
        self.assertEqual(expected, result)

    def testSumRates(self):
        expected = 0.004637967436947257
        result = self.sim_on_branch1.sum_rates()
        self.assertEqual(expected, result)

        expected = 0.00954128577225131
        result = self.sim_on_branch2.sum_rates()
        self.assertEqual(expected, result)

        expected = 0.013399271035090315
        result = self.sim_on_branch3.sum_rates()
        self.assertEqual(expected, result)

        expected = 0.0030337488758443035
        result = self.sim_on_branch4.sum_rates()
        self.assertEqual(expected, result)


# ==========================================
# Tests for SimulateOnTree
# ==========================================
class TestSimulateOnTree(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        random.seed(9001)  # Set seed for pseudo-random number generator

        s1 = 'GTACGATCGATCGATGCTAGC'
        kappa = 0.3
        mu = 0.0005
        pi1 = Sequence.get_frequency_rates(s1)
        omegas = [0.29327471612351436, 0.6550136761581515, 1.0699896623909886, 1.9817219453273531]
        sequence1 = Sequence(s1, {'+0': [(0, 20)]}, kappa, mu, pi1, omegas)

        phylo_tree = Phylo.read(TEST_TREE, 'newick', rooted=True)
        self.sim_on_tree1 = SimulateOnTree(sequence1, phylo_tree)

    def testGetParentClade(self):
        # Label each clade
        ch = 'a'
        for child_clade in self.sim_on_tree1.phylo_tree.find_clades(order='level'):
            child_clade.name = ch
            ch = chr(ord(ch) + 1)

        results = []
        for child_clade in self.sim_on_tree1.phylo_tree.find_clades(order='level'):
            result = self.sim_on_tree1.get_parent_clade(child_clade)
            print(result)
            results.append(result.name)

        expected = ['a', 'a', 'a', 'a', 'd', 'd']
        self.assertEqual(expected, results)


if __name__ == '__main__':
    unittest.main()
