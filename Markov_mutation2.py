import random

# # Compute the nucleotide frequencies in a DNA string
def get_base_frequencies_v2(dna):
    return {base: dna.count(base) / float(len(dna))
            for base in 'ATGC'}

def format_frequencies(frequencies):
    return ', '.join(['%s: %.2f' % (base, frequencies[base])
                      for base in frequencies])

# Draw a random letter according to the probabilities from markov_chain[b], where v us tge vase at current position
def draw(discrete_probdist):
    """
    Draw random value from discrete probability distribution
    represented as a dict: P(x=value) = discrete_probdist[value].
    """
    # Method:
    # http://en.wikipedia.org/wiki/Pseudo-random_number_sampling
    limit = 0
    r = random.random()
    for value in discrete_probdist:
        limit += discrete_probdist[value]
        if r < limit:
            return value


# Create a transition probability matrix
def create_markov_chain():
    markov_chain = {}
    for from_base in 'ATGC':
        # Generate random transition probabilities by dividing
        # [0,1] into four intervals of random length
       slice_points = sorted(
           [0] + [random.random()for i in range(3)] + [1])
       transition_probabilities = \
           [slice_points[i+1] - slice_points[i] for i in range(4)]
       markov_chain[from_base] = {base: p for base, p
                         in zip('ATGC', transition_probabilities)}
    return markov_chain

mc = create_markov_chain()
#print (mc)
#print (mc['A']['T']) # probability of transition from A to T

# Mutate a nucleotide in a random site to a base according to the probability matrix
def mutate_via_markov_chain(dna, markov_chain):
    dna_list = list(dna)
    mutation_site = random.randint(0, len(dna_list) - 1)
    from_base = dna[mutation_site]
    to_base = draw(markov_chain[from_base])
    dna_list[mutation_site] = to_base
    return ''.join(dna_list)



# Testing code
dna = 'TTACGGAGATTTCGGTATGCATATGGTGCCATGA'
print ('Starting DNA:', dna)
print (format_frequencies(get_base_frequencies_v2(dna)))

nmutations = 10000
for i in range(nmutations):
    dna = mutate_via_markov_chain(dna, mc)

print ('DNA after %d mutations (Markov chain):' % nmutations, dna)
print (format_frequencies(get_base_frequencies_v2(dna)))