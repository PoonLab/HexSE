from evol_rates import get_reading_frames
from evol_rates import codon_iterator
#from evol_rates import get_omega
from evol_rates import draw_omega_values
from evol_rates import get_codon
from evol_rates import get_syn_subs
from evol_rates import Rates
from evol_rates import reverse_and_complement
from evol_rates import update_rates
from sequence_info import sort_orfs, Nucleotide, Sequence
from evol_rates import update_rates



bias = { 'A': {'C':0.001, 'G':0.065 ,'T':0.002 },
         'C': {'A':0.001, 'G':0.00001 ,'T':0.064 },
         'G': {'A':0.065, 'C':0.00001 ,'T':0.00001 },
         'T': {'A':0.002, 'C':0.064 ,'G':0.00001 } }

mu = 1
pi = None
omega = None

#
original_seq = 'AAAAAAAGATTTAAAA'
orfs = [(2,10), (4,15), (8,0), (10,2), (9,3)] #, (0,5)
sorted_orfs = sort_orfs(orfs)
seq = Sequence(original_seq, sorted_orfs)
print(seq.codon)
#rates = Rates(seq, mu, sorted_orfs, bias, pi, omega)
#print(rates.before_omega)
#print(rates.seq)

#n = update_rates(rates, 5, "A")

print(get_codon(original_seq, 5, orfs[1]))