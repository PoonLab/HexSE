# Calculate evolution rates across the sequence

transition_score = -5
transversion_score = -1
synonymous_score = 5
nonsynonymous_score = 1
rate_matrix = {
    'A': {
        'A': 0,
        'C': transversion_score,
        'T': transversion_score,
        'G': transition_score,
    },

}

def get_frequency_rates(seq):
    seq = list(seq)
    frequencies = {
        'A': 0,
        'C': 0,
        'T': 0,
        'G': 0,
    }

    for nucleotide in frequencies:
        frequencies[nucleotide] = seq.count(nucleotide)

    return { nucleotide: ocurrences/len(seq) for nucleotide, ocurrences in frequencies.items() }


def get_syn_nonsyn_score():
    return 1


def get_tt_score(current, mutated):
    return 1

# bias = { 'A': [C,G,T],
#          'C': [A,G,T],
#          'G': [A,C,T],
#          'T': [A,C,G],
# }


def evol_rates(seq, mu, pi, bias, omega):
    """
    Generate rate vector from sequence given parameters.
    @param seq: the nucleotide sequence of length <L>.  Assumed to start and end in reading frame +0.
    @param mu: the global rate (substitutions/site/unit time).
    @param pi: <optional> a vector of stationary nucleotide frequencies.  If not specified, defaults
               to the empirical frequencies in <seq>.
    @param bias: <optional> a List of 6 rates [AC, AG, AT, CG, CT, GT] assuming time-reversibility.
                    If not specified, defaults to 1's.
    @param omega: <optional> a List of dN/dS (omega) ratios along sequence length as {dict}.
                  {dict} always contains key +0 (parent reading frame).  May contain omega for alternate
                  reading frame as (+1, +2, -0, -1 or -2).  Codon position is determined by the nt's
                  location relative to start of <seq>.
    @return rates: a List of tuples for rates of 3 possible nucleotide substitution at each site
                   in <seq> (3xL).
    """
    L = len(seq)
    rates = [(mu, mu, mu) for i in range(L)]  # initialize our list with baseline rate
    frequency_rates = get_frequency_rates(seq)
    # iterate over every nt in seq
    for i in range(L):
    # 1. what is the current nucleotide X?
    # 2. apply stationary frequencies (pi) to tuple given X
        pi = frequency_rates(seq[i])
        rates[i] = tuple([pi*j for j in rates[i]])
    # 3. apply biases to tuple
        if seq[i] == 'A':
            for j in rates[i]:
                rates[i][j] = rates[i][j] * bias['A'][j]

        elif seq[i] == 'C':
            for j in rates[i]:
                rates[i][j] = rates[i][j] * bias['C'][j]

        elif seq[i] == 'G':
            for j in rates[i]:
                rates[i][j] = rates[i][j] * bias['G'][j]

        elif seq[i] == 'T':
            for j in rates[i]:
                rates[i][j] = rates[i][j] * bias['T'][j]

    # 4. apply omega in parent reading frame
    # 5. if alternate reading frame(s) is present, apply other omega(s)
    
    # Step 1: Calculate frequency rate of nucleotides
    seq = list(seq)
    # Step 2: Calculate mutation rate for each nucleotide in the sequence
    for nucleotide in seq:
        mutation_rates.append({
            'A': frequency_rates[nucleotide] * get_tt_score(nucleotide, 'A') * get_syn_nonsyn_score(),
            'C': frequency_rates[nucleotide] * get_tt_score(nucleotide, 'C') * get_syn_nonsyn_score(),
            'G': frequency_rates[nucleotide] * get_tt_score(nucleotide, 'G') * get_syn_nonsyn_score(),
            'T': frequency_rates[nucleotide] * get_tt_score(nucleotide, 'T') * get_syn_nonsyn_score(),
        })

    # Map of nucleotides with their respective mutation rate
    return mutation_rates

# if __name__ == '__main__':
#     print(evol_rates(('AA')))
