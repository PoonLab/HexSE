from evol_rates import evol_rates
from evol_rates import get_frequency_rates

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

print(evol_rates(seq, mu, bias, pi, omega))
#seq_rates = evol_rates(seq, mu, bias, pi = None)
