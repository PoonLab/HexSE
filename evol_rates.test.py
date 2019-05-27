from evol_rates import get_reading_frames
from evol_rates import codon_iterator
#from evol_rates import get_omega
from evol_rates import get_codon
from evol_rates import get_syn_subs
from evol_rates import Rates
from evol_rates import reverse_and_complement
from evol_rates import update_rates
from sequence_info import sort_orfs, Nucleotide, Sequence


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
original_seq = 'AAAAAAGAAAAA'
orfs = [(2,10), (4,15), (0,8), (9,4), (8,0), (10,2)]
sorted_orfs = sort_orfs(orfs)
# seq = []
# for position in range(len(seq_original)):
#     seq.append(Nucleotide(seq_original[position], position))
# seq = Sequence(original_seq, sorted_orfs)
# print(seq)
# print(seq[3].orf_presence)


#rates = Rates(seq, mu, bias, pi, orfs)
#(update_rates(rates, 6, 'C'))
#print(rates.pi)
#print(get_omega(orfs))


# print(sorted_orfs)
# n = Nucleotide('A', 3, sorted_orfs)
# print(n.orf_presence)
#
#
#update_rates [{}]
# Class Nuclotide: position, letra, tupla (TRUE, TRUE, FALSE, FALSE, FALSE)
