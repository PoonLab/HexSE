from ovrf_functions import reverse_and_complement
from ovrf_functions import sort_orfs
from sequence_info import Sequence
from evol_rates import Rates
from evol_rates import update_rates
from simulate import Simulate
from Bio import Phylo

bias = { 'A': {'C':0.001, 'G':0.065 ,'T':0.002 },
         'C': {'A':0.001, 'G':0.00001 ,'T':0.064 },
         'G': {'A':0.065, 'C':0.00001 ,'T':0.00001 },
         'T': {'A':0.002, 'C':0.064 ,'G':0.00001 } }

mu = 1
pi = None
omega = None

# Inputs
original_seq = 'AAAAAAAAAAAAAAATTTTA'
orfs = [(2,10), (4,15), (8,0), (10,2)] #, (0,5), (9,3)


sorted_orfs = sort_orfs(orfs)
seq = Sequence(original_seq, sorted_orfs)
rates = Rates(seq, mu, bias, pi, omega)

n = update_rates(rates, 8, "T")

s = Simulate(rates)
tree = Phylo.read("test_tree.txt", "newick", rooted = True)
print(s.traverse_tree(tree))