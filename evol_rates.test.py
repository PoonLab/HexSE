from evol_rates import get_reading_frames
from evol_rates import codon_iterator
from evol_rates import get_omega
from evol_rates import get_codon
from evol_rates import get_syn_subs
from evol_rates import Rates
import numpy as np
from evol_rates import reverse_and_complement
from evol_rates import update_rates


bias = { 'A': {'C':0.001, 'G':0.065 ,'T':0.002 },
         'C': {'A':0.001, 'G':0.00001 ,'T':0.064 },
         'G': {'A':0.065, 'C':0.00001 ,'T':0.00001 },
         'T': {'A':0.002, 'C':0.064 ,'G':0.00001 } }

mu = 1
seq = 'AGTCGTGCTTCGG'
pi = None
orfs = ((0,8),(1,9))

codon_dict = {  'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
                '---':'-', 'XXX':'?'}
#
seq = 'AAAAAAAAAA'
#orfs = ((2,10),(10,2))

#print(get_syn_subs(seq,orfs,codon_dict), "\n")
rates = Rates(seq, mu, bias, pi, orfs)
# print(rates)
print(update_rates(rates, seq, 6, 'C', (2,10)))


