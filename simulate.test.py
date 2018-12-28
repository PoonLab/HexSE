import unittest
from unittest.mock import patch
from simulate import Simulate

class TestSimulate(unittest.TestCase):
    def test_simulate_on_branch_generic(self):
        sim = Simulate(rate=None, matrix=None)
        alphabet = list(sim.alphabet)
        mutated = sim.simulate_on_branch('AAA', 0.3)
        self.assertEqual(len(mutated), 3)
        for nucleotide in mutated:
            self.assertIn(nucleotide, alphabet)

    @patch('simulate.random')
    def test_simulate_on_branch_one_mutation(self, random):
        # Mock random values to force a resulting mutated sequence
        # random 0.3 to mutate just once
        # random > 0.075 to select A.C from the matrix (0.825)
        random.random.return_value = 0.3
        # randint 2 to select the third position
        random.randint.return_value = 2

        matrix = {
            'A': {
                'A': 0.025,
                'T': 0.025,
                'G': 0.025,
                'C': 0.825,
            },
            # The rest of the matrix nucleotides...
        }

        sim = Simulate(None, matrix)
        alphabet = list(sim.alphabet)
        mutated = sim.simulate_on_branch('AAAA', 0.3)
        self.assertEqual(mutated, 'AACA')

    @patch('simulate.random')
    def test_simulate_on_branch_multiple_mutations(self, random):
        # random 0.1 to mutate three times
        # random > 0.075 to select A.C from the matrix (0.825)
        random.random.return_value = 0.1
        # randint returns 3, 4, 5 in each mutation to get the last
        # three nucleotides mutated
        random.randint.side_effect = [3, 4, 5]

        matrix = {
            'A': {
                'A': 0.025,
                'T': 0.025,
                'G': 0.025,
                'C': 0.825,
            },
            # The rest of the matrix nucleotides...
        }

        sim = Simulate(None, matrix)
        alphabet = list(sim.alphabet)
        mutated = sim.simulate_on_branch('AAAAAA', 0.3)
        self.assertEqual(mutated, 'AAACCC')

if __name__ == '__main__':
    unittest.main()
