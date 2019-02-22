from evol_rates import evol_rates
from evol_rates import get_frequency_rates

# Substitution bias (transition-transversion) AC:0.001, AG:0.065, AT:0.002, CG:0.00001, CT:0.064, GT:0.00001
# Values for this matrix were taken from: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0070388

bias = { 'A': [0.001,0.065,0.002],
         'C': [0.001,0.00001,0.064],
         'G': [0.065,0.00001,0.00001],
         'T': [0.002,0.064,0.00001], }


mu = 1

seq = 'AACTGAA'

#print(evol_rates(seq, mu, bias,pi =  None))
seq_rates = evol_rates(seq, mu, bias, pi=None)


# Hypothetical dictionary constructed for an 'AAA' (lys) sequence, where every mutation in first and second position are nonsyn
omega = {'reading_frame1': [{'C': 0.001, 'G': 0.001, 'T':0.001}, {'C': 0.001, 'G': 0.001, 'T':0.001},{'C': 0.001, 'G': 0.3, 'T':0.001}],
         'reading_frame2': [{'C': 0.001, 'G': 0.001, 'T':0.001}, {'C': 0.001, 'G': 0.001, 'T':0.001},{'C': 0.001, 'G': 0.3, 'T':0.001}]}
print(omega)

for key, value in omega.items():
    second_ratio_dictionary = value[1]
    c_ratio_number = second_ratio_dictionary['C']
    print(c_ratio_number)
    # print(value[1]['C'])