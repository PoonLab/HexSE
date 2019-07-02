from ovrf_functions import reverse_and_complement
from ovrf_functions import sort_orfs
from sequence_info import Sequence
from evol_rates import Rates
from simulate import Simulate
from Bio import Phylo, SeqIO
from orfs_from_genbank import orfs_from_gb
from orfs_from_genbank import sq_from_gb

bias = { 'A': {'C':0.001, 'G':0.065 ,'T':0.002 },
         'C': {'A':0.001, 'G':0.00001 ,'T':0.064 },
         'G': {'A':0.065, 'C':0.00001 ,'T':0.00001 },
         'T': {'A':0.002, 'C':0.064 ,'G':0.00001 } }

mu = 0.00001
pi = None
omega = None

gb_file = "NC_003977.2_HBV.gb"
# Inputs
original_seq = sq_from_gb(gb_file)
orfs = orfs_from_gb(gb_file)
print(orfs)
#orfs = [(2,10), (4,9), (8,0), (10,2)] #, (0,5), (9,3)


sorted_orfs = sort_orfs(orfs)
print(sorted_orfs)

seq = Sequence(original_seq, sorted_orfs)
rates = Rates(seq, mu, bias, pi, omega)
#print(rates.mu)
tree = Phylo.read("/home/lmunoz/Projects/ovrf/HBV/abayomi_tree.nwk", "newick", rooted = True)
s = Simulate(rates, tree)
print(s.sum_rates())

print(s.get_alignment())