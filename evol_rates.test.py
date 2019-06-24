from ovrf_functions import reverse_and_complement
from ovrf_functions import sort_orfs
from sequence_info import Sequence
from evol_rates import Rates
from simulate import Simulate
from Bio import Phylo, SeqIO

bias = { 'A': {'C':0.001, 'G':0.065 ,'T':0.002 },
         'C': {'A':0.001, 'G':0.00001 ,'T':0.064 },
         'G': {'A':0.065, 'C':0.00001 ,'T':0.00001 },
         'T': {'A':0.002, 'C':0.064 ,'G':0.00001 } }

mu = 1
pi = None
omega = None

# Inputs
record = SeqIO.read("HBV.fasta", "fasta")
original_seq = str(record.seq)
orfs = [(2,10), (4,9), (8,0), (10,2)] #, (0,5), (9,3)


sorted_orfs = sort_orfs(orfs)
seq = Sequence(original_seq, sorted_orfs)
rates = Rates(seq, mu, bias, pi, omega)
tree = Phylo.read("abayomi_tree.nwk", "newick", rooted = True)
s = Simulate(rates, tree)
print(s.get_alignment())
