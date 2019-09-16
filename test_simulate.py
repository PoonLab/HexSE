from refact.new_sequence_info import *
import unittest
from refact.simulation_refactore import Simulate
from Bio import Phylo


class testSumRates(unittest.TestCase):

    def testSimpleRates(self):

        s = Sequence('ATGAAAGGGTTTACATAG')
        mu = 0.00001
        bias = {'A': {'C': 0.001, 'G': 0.065, 'T': 0.002},
                'C': {'A': 0.001, 'G': 0.00001, 'T': 0.064},
                'G': {'A': 0.065, 'C': 0.00001, 'T': 0.00001},
                'T': {'A': 0.002, 'C': 0.064, 'G': 0.00001}}
        pi = None
        omega = None
        r = Rates(s, mu, bias, pi, omega)
        tree = '(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);'
        sim = Simulate(r, tree)
        result = sim.sum_rates()
        expected = 3.7298260000000006e-06
        self.assertEqual(expected, result)