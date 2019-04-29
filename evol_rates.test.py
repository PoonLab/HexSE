from evol_rates import evol_rates
from evol_rates import get_reading_frames
from evol_rates import codon_iterator
from evol_rates import get_omega
from evol_rates import get_codon
from evol_rates import get_syn_subs

# Substitution bias (transition-transversion) AC:0.001, AG:0.065, AT:0.002, CG:0.00001, CT:0.064, GT:0.00001
# Values for this matrix were taken from: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0070388

# bias = { 'A': [C,G,T],
#          'C': [A,G,T],
#          'G': [A,C,T],
#          'T': [A,C,G],
# }

bias = { 'A': [0.001,0.065,0.002],
         'C': [0.001,0.00001,0.064],
         'G': [0.065,0.00001,0.00001],
         'T': [0.002,0.064,0.00001], }
mu = 1
seq = 'AAA'
pi = None

# Hypothetical dictionary constructed for an 'AAA' (lys) sequence, where every mutation in first and second position are nonsyn
# omega = {'reading_frame0plus': [[ 0.001,  0.001, 0.001], [ 0.001,  0.001, 0.001],[ 0.001,  0.3, 0.001]],
#          'reading_frame1plus': [[ 0.001,  0.001, 0.001], [ 0.001,  0.001, 0.001],[ 0.001,  0.3, 0.001]]}

# draw omega values from some probability distribution (e.g., gamma)
#  or user-specified vectors
# omega = {'+0': [0.1, 0.13, 0.09, 1.5, 0.7, 0.12],
#          '+1': [2.7, 0.2, 0.17, 0.56, 0.2, 0.3]}

#print(evol_rates(seq, mu, bias, pi, omega))
#seq_rates = evol_rates(seq, mu, bias, pi = None)

#Test get_reading_frames
#seq='ATG*ATG*ATG***TAG*TAG*TAG'
seq='AT GTA TGT ATG TTTTAGTTAGTTAG'
seq = 'AAAAAATTTTTAAAAGGGATATAGATAC'

#Provided by the user
#reading_frames = get_reading_frames(seq)
# print(reading_frames)
#
#
#print(get_omega(reading_frames))
#
# omega = get_omega(reading_frames)0
# print(type(omega))
# print(evol_rates(seq, mu, bias, pi, reading_frames))
# print(get_syn_codons('ATT'))

# Testing get_codos
#orf = [2,13]
#info = get_codon(seq, 8, orf)
#print(info[0])
#print(info[1])

# Testing get_syn_subs:

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

seq = 'AAAAAAAAAAAA'
orfs = [[0,8],[1,9],[2,7]]
#print(len(orfs))
orfs = ((0,8),(1,9))
print("IS IT SIN:", get_syn_subs(seq,orfs,codon_dict), "\n")
print("EVOL RATES:", evol_rates(seq, mu, bias, pi, orf), "\n")
print("OMEGA:", (get_omega(orfs)), "\n")
