from evol_rates import evol_rates
from evol_rates import get_reading_frames


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
omega = {'reading_frame0plus': [[ 0.001,  0.001, 0.001], [ 0.001,  0.001, 0.001],[ 0.001,  0.3, 0.001]],
         'reading_frame1plus': [[ 0.001,  0.001, 0.001], [ 0.001,  0.001, 0.001],[ 0.001,  0.3, 0.001]]}

# draw omega values from some probability distribution (e.g., gamma)
#  or user-specified vectors
omega = {'+0': [0.1, 0.13, 0.09, 1.5, 0.7, 0.12], 
         '+1': [2.7, 0.2, 0.17, 0.56, 0.2, 0.3]}

#print(evol_rates(seq, mu, bias, pi, omega))
#seq_rates = evol_rates(seq, mu, bias, pi = None)

#Test get_reading_frames
seq='AAATGBBBATGTAG***TAG**'
print (len(seq))
print(get_reading_frames(seq))




