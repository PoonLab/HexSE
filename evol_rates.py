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
        frequencies[nucleotide] = round((float(seq.count(nucleotide))/(len(seq))),2)


    return frequencies


def get_syn_nonsyn_score():
    return 1


def get_tt_score(current, mutated):
    return 1

# bias = { 'A': [C,G,T],
#          'C': [A,G,T],
#          'G': [A,C,T],
#          'T': [A,C,G],
# }

def evol_rates(seq, mu, bias, pi):
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
    @return seq_rates: a List of tuples for rates of 3 possible nucleotide substitution at each site
                   in <seq> (3xL).
    """

    L = len(seq)
    seq_rates = [(mu, mu, mu) for position in range(L)]  # initialize our list with baseline rate
    frequency_rates = get_frequency_rates(seq)
    # iterate over every nt in seq
    for position in range(L):
        # 1. what is the current nucleotide X?
        nucleotide = seq[position]

        # 2. apply stationary frequencies (pi) to tuple given X
        if pi is None:
            pi = frequency_rates[nucleotide]

        seq_rates[position] = tuple([pi*j for j in seq_rates[position]])
        # 3. apply biases to tuple
        biased_rates = []

        for rate in range(3):
            local_rates = seq_rates[position][rate] * bias[nucleotide][rate]
            biased_rates.append(local_rates)

        formatted_biased_rates = [ '%.6f' % elem for elem in biased_rates]
        seq_rates[position] = tuple(formatted_biased_rates)

    # 4. apply omega in parent reading frame



    # 5. if alternate reading frame(s) is present, apply other omega(s)


    return seq_rates
