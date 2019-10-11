# Store sequence information

import scipy
import scipy.stats as ss
import numpy as np
import random

from event_tree import EventTree

TRANSITIONS_DICT = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}

NUCLEOTIDES = ['A', 'C', 'G', 'T']

COMPLEMENT_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'W', 'R': 'Y', 'K': 'M', 'Y': 'R',
                   'S': 'S', 'M': 'K', 'B': 'V', 'D': 'H',
                   'H': 'D', 'V': 'B', '*': '*', 'N': 'N',
                   '-': '-'}

CODON_DICT = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
              'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
              'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
              'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
              'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
              'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
              'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
              '---': '-', 'XXX': '?'}


class Sequence:
    """
    Store inputs and create sequence objects
    :ivar original_seq: the nucleotide sequence as a string
    :ivar nt_sequence: list of Nucleotide objects
    """

    def __init__(self, original_seq, rcseq, sorted_orfs, mu, pi, kappa, omega=None):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.
        :param original_seq: The nucleotide sequence as a string
        :param rcseq: The reverse complement of the sequence
        :param sorted_orfs: A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        :param mu: The global rate (substitutions/site/unit time)
        :param kappa: transition transversion rate ratio
        :param pi: a vector of stationary nucleotide frequencies.
        :param omega: <optional> a List of dN/dS (omega) ratios.
                     If user does not specify, they are drawn from a gamma distribution
        """

        self.original_seq = original_seq
        self.orfs = {}
        self.rcseq = rcseq
        self.orfs = sorted_orfs
        self.mu = mu
        self.pi = pi
        self.kappa = kappa
        self.omega_values = omega

        # Drawn omega values from gamma distribution
        if self.omega_values is None:
            self.omega_values = self.get_omega_values(2, 4)

        # Calculate stationary frequency rates
        if self.pi is None:
            self.pi = self.get_frequency_rates(self.original_seq)

        # Create event tree for sequence (tree containing all possible mutation and the parameters that should be taken into account when calculating rates)
        self.tree = EventTree(self.pi, self.kappa)
        self.event_tree = self.tree.create_event_tree()

        # Create Nucleotides
        self.nt_sequence = DoubleLinkedList()
        for pos_in_seq, nt in enumerate(self.original_seq):
            self.nt_sequence.insert_nt(nt, pos_in_seq)

        # Set Codons based on the reading frames
        if self.orfs is not None:
            for frame in self.orfs:
                orf_list = self.orfs[frame]
                for orf in orf_list:
                    codons = self.find_codons(frame, orf)

                    # Tell Nucleotide which Codon(s) it belongs to
                    for codon in codons:
                        for i, nt in enumerate(codon.nts_in_codon):
                            nt.codons.append(codon)

        # Calculate mutation rates for each nucleotide in sequence, populate the event tree which each nucleotide
        for nt in iter(self.nt_sequence):
            rates = self.get_substitution_rates(nt)
            nt.set_rates(rates[0])
            nt.set_my_omegas(rates[1])
            nt.set_nt_rate()

        # Update event_tree to include a list of nucleotides on the tips
        self.event_tree = self.get_nts_on_tips()

    def get_nts_on_tips(self):
        """
        Look at the tips of the event tree and create keys and values for nucleotides in substitution and number of events
        :param event_tree: Dictionary of class EventTree populated with all of the nucleotides associated with each event
        :return: Updated event Tree that includes: - List of nucleotides for each tip
                                                   - Total number of events
                                                   - Number of events leading to each to_nt
        """
        updated_event_tree = self.event_tree
        my_tree = updated_event_tree['to_nt']

        total_events = 0  # Number of all possible events on the tree
        for key1, to_nt in my_tree.items():
            subset = to_nt['from_nt']
            events_for_to_nt = 0  # Number of all events that can lead to a mutation
            for key2, from_nt in subset.items():
                if from_nt:
                    nt_in_substitution = []  # Nucleotides associated with each substitution event
                    # Add nucleotides that are not involved in any non-syn substitution
                    if from_nt['is_syn']:
                        nt_in_substitution.extend([nts for nts in from_nt['is_syn']])
                    # Add nucleotides involved in non-syn substitutions
                    non_syn_subs = from_nt['is_nonsyn']
                    for key3, nts in non_syn_subs.items():
                        nt_in_substitution.extend([nt for nt in nts])
                    updated_event_tree['to_nt'][key1]['from_nt'][key2].update([('nts_in_subs', nt_in_substitution)])
                    updated_event_tree['to_nt'][key1]['from_nt'][key2].update([('number_of_events', len(nt_in_substitution))])
                    events_for_to_nt += len(nt_in_substitution)
                    total_events += len(nt_in_substitution)

            updated_event_tree['to_nt'][key1].update([('events_for_nt', events_for_to_nt)])

        updated_event_tree['total_events'] = total_events
        return updated_event_tree

    def get_substitution_rates(self, nt):
        """
        Calculates substitution rates for each mutation and populates Event tree
        :param nt: object of class Nucleotide
        :return: 1. Dictionary of substitutions rates, keyed by nt subs
                 2. Dictionary with omega keys
        """
        sub_rates = {}
        my_omega_keys = {}
        current_nt = nt.get_state()

        for to_nt in NUCLEOTIDES:
            if to_nt == current_nt:
                sub_rates[to_nt] = None
                my_omega_keys[to_nt] = None
            else:
                # Apply global substitution rate and stationary nucleotide frequency
                sub_rates[to_nt] = self.mu * self.pi[current_nt]
                if self.is_transv(current_nt, to_nt):
                    sub_rates[to_nt] *= self.kappa

                chosen_omegas = [0, 0, 0, 0]      # omegas applied given a substitution from current_nt to to_nt
                for codon in nt.codons:
                    pos_in_codon = codon.nt_in_pos(nt)
                    if codon.is_nonsyn(pos_in_codon, to_nt):    # Apply omega when mutation is non-synonymous
                        omega_index = random.randrange(len(self.omega_values))
                        sub_rates[to_nt] *= self.omega_values[omega_index]
                        chosen_omegas[omega_index] += 1

                omegas_in_sub = tuple(chosen_omegas)
                # key = self.create_keys(omegas_in_subs)
                my_omega_keys[to_nt] = omegas_in_sub
                # Populate even tree using omega_keys
                if any(omegas_in_sub):  # At least one omega was used
                    current_event = self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']
                    # Create key if needed, associate it with the current nucleotide
                    if omegas_in_sub not in current_event:
                        self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn'][omegas_in_sub] = [nt]
                    else:
                        self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn'][omegas_in_sub].append(nt)

                else:  # If mutation is syn in all codons
                    self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_syn'].append(nt)

        return sub_rates, my_omega_keys

    @staticmethod
    def is_transv(from_nt, to_nt):

        if from_nt == to_nt:
            transv = None
        else:
            transv = True
            if TRANSITIONS_DICT[from_nt] == to_nt:
                transv = False
        return transv

    def get_sequence(self):
        return self.nt_sequence

    @staticmethod
    def codon_iterator(my_orf, start_pos, end_pos):
        """
        Generator to move every three nucleotides (codon)
        :param my_orf: A list of Nucleotides in the ORF
        :param start_pos: The start position of the ORF
        :param end_pos: The end position of the ORF
        :yield codon
        """
        if start_pos > end_pos:  # Negative strand
            my_orf.reverse()

        i = 0
        while i < len(my_orf):
            yield my_orf[i:i+3]
            i += 3

    def find_codons(self, frame, orf):
        """
        Gets the Codon sequence
        :param frame: the frame of the ORF
        :param orf: tuple containing the coordinates of the ORF
        :return: a list of Codon objects for the specified ORF
        """
        codons = []
        start_pos = orf[0]
        end_pos = orf[1]
        my_orf = self.nt_sequence.slice_sequence(start_pos, end_pos)

        # Iterate over list by threes and create Codons
        for cdn in self.codon_iterator(my_orf, start_pos, end_pos):
            codon = Codon(frame, orf, cdn)
            codons.append(codon)

        return codons

    @staticmethod
    def get_frequency_rates(seq):
        """
        Frequency of nucleotides in the DNA sequence
        :param seq: the DNA sequence
        :return: a dictionary frequencies where the key is the nucleotide and the value is the frequency
        """
        frequencies = {
            'A': 0,
            'C': 0,
            'T': 0,
            'G': 0
        }

        for nucleotide in frequencies:
            frequencies[nucleotide] = round((float(seq.count(nucleotide)) / (len(seq))), 2)

        return frequencies

    @staticmethod
    def get_omega_values(alpha, ncat):
        """
        Draw ncat number of omega values from a discretized gamma distribution
        :param alpha: shape parameter
        :param ncat: Number of categories (expected omegas)
        :return: list of ncat number of omega values (e.i. if ncat = 3, omega_values = [0.29, 0.65, 1.06])
        """
        values = Sequence.discretize_gamma(alpha=alpha, ncat=ncat)
        omega_values = list(values)
        return omega_values

    @staticmethod
    def discretize_gamma(alpha, ncat, dist=ss.gamma):
        """
        Divide the gamma distribution into a number of intervals with equal probability and get the mid point of those intervals
        From https://gist.github.com/kgori/95f604131ce92ec15f4338635a86dfb9
        :param alpha: shape parameter
        :param ncat: Number of categories
        :param dist: function from scipy stats
        :return: array with ncat number of values
        """
        if dist == ss.gamma:
            dist = dist(alpha, scale=1 / alpha)
        elif dist == ss.lognorm:
            dist = dist(s=alpha, scale=np.exp(0.5 * alpha ** 2))
        quantiles = dist.ppf(np.arange(0, ncat) / ncat)
        rates = np.zeros(ncat, dtype=np.double)  # return a new array of shape ncat and type double
        for i in range(ncat - 1):
            rates[i] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x),
                                                   quantiles[i], quantiles[i + 1])[0]
        rates[ncat - 1] = ncat * scipy.integrate.quad(lambda x: x * dist.pdf(x),
                                                      quantiles[ncat - 1], np.inf)[0]
        return rates

    def get_string_sequence(self):
        seq_as_string = ''
        for nt in iter(self.nt_sequence):
            seq_as_string += nt.state
        return seq_as_string


class Nucleotide:
    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    :ivar state: the nucleotide base (A, T, G, C)
    :ivar pos_in_seq: the position of the nucleotide relative to the start of the sequence
    :ivar left_nt: reference to the nucleotide's left neighbour
    :ivar right_nt: reference to the nucleotide's right neighbour
    """

    def __init__(self, state, pos_in_seq, left_nt=None, right_nt=None):

        """
        :param state: Nucleotide A, C, G or T
        :param pos_in_seq: Position of the nucleotide in the sequence
        :param left_nt : reference to the adjacent nucleotide to the left (default to None)
        :param right_nt: reference to the adjacent nucleotide to the right (default to None)
        """
        self.state = state
        self.pos_in_seq = pos_in_seq
        self.left_nt = left_nt
        self.right_nt = right_nt
        self.codons = []
        self.complement_state = COMPLEMENT_DICT[self.state]
        self.rates = {}  # Mutation rates
        self.my_omegas = []  # omegas that have used for this nucleotide when creating rates

    def __repr__(self):
        return self.state

    def get_state(self):
        return self.state

    def get_pos_in_seq(self):
        return self.pos_in_seq

    def get_left_nt(self):
        return self.left_nt

    def get_right_nt(self):
        return self.right_nt

    def get_complement_state(self):
        return self.complement_state

    def set_state(self, new_state):
        self.state = new_state

    def set_pos_in_seq(self, new_pos):
        self.pos_in_seq = new_pos

    def set_left_nt(self, new_left_nt):
        self.left_nt = new_left_nt

    def set_right_nt(self, new_right_nt):
        self.right_nt = new_right_nt

    def add_codon(self, codon):
        self.codons.append(codon)

    def set_rates(self, rates):
        self.rates = rates

    def set_my_omegas(self, omegas):
        self.my_omegas = omegas

    def get_mutation_rate(self):
        total_rate = 0
        for to_nt, value in self.rates.items():
            if value:
                total_rate += value
        return total_rate

    def set_complement_state(self):
        self.complement_state = COMPLEMENT_DICT[self.state]

    def set_nt_rate(self):
        self.mutation_rate = self.get_mutation_rate()


class DoubleLinkedList:
    """
    Double linked list linking together objects of class Nucleotide
    default initialization with empty head node
    """

    def __init__(self):
        self.head = None  # head node (starting nucleotide)
        self.current_nt = None  # Pointer to current nt for insertion
        self.next_iter_nt = None  # Current state of the iteration

    def __iter__(self):
        self.next_iter_nt = self.head
        return self

    def __next__(self):
        # Note: Not thread safe
        if self.next_iter_nt is not None:
            nt = self.next_iter_nt
            self.next_iter_nt = nt.get_right_nt()
            return nt
        else:
            raise StopIteration

    def insert_nt(self, state, position):
        """
        Insert objects of class Nucleotide to the end of the DoubleLinkedList
        :param state: Nucleotide state in sequence
        :param position: Position of nt in sequence
        """
        new_nt = Nucleotide(state, position)  # create new Nucleotide object

        # Assign the first nucleotide as head
        if self.head is None:
            self.head = new_nt
            self.current_nt = new_nt

        else:
            new_nt.set_left_nt(self.current_nt)  # For the new nucleotide, create a left pointer towards the current one
            self.current_nt.set_right_nt(new_nt)  # Create the double link between current and new
            self.current_nt = new_nt

    def get_head(self):
        return self.head

    def print_seq(self):  # Print the string of nucleotides (check the class is working properly)
        s = ''
        temp = self.head

        while temp is not None:
            s += temp.get_state()
            temp = temp.get_right_nt()
        print(s)

    def nucleotide_at_pos(self, position):
        """
        Traverse sequence to find the Nucleotide object at a specific position
        :param position: the position of the Nucleotide
        :return: the Nucleotide object in the specified position
        """
        current_nt = self.get_head()
        while current_nt is not None:
            if position == current_nt.get_pos_in_seq():
                return current_nt
            current_nt = current_nt.get_right_nt()
        return None

    def slice_sequence(self, start_pos, end_pos):
        """
        Slices the Nucleotide sequence from the start position, up to and including the end position (inclusive slicing)
        :param start_pos: the start position
        :param end_pos: the end position
        :return sub_seq: a list of Nucleotides between the start and end positions (in 5', 3' direction)
        """

        sub_seq = []
        if start_pos < end_pos:  # Positive strand
            curr_nt = self.nucleotide_at_pos(start_pos)
        else:   # Negative strand
            curr_nt = self.nucleotide_at_pos(end_pos)
            end_pos = start_pos

        while curr_nt is not None and curr_nt.pos_in_seq <= end_pos:
            sub_seq.append(curr_nt)
            curr_nt = curr_nt.right_nt

        return sub_seq


class Codon:
    """
    Stores information about the frameshift, ORF, omega, and pointers to 3 Nucleotide objects
    """

    def __init__(self, frame, orf, ptrs_to_nts):
        """
        Create a Codon
        :param frame: the reading frame (+0, +1, +2, -0, -1, -2)
        :param orf: a tuple containing the reading frame and the coordinates of the orf
        :param ptrs_to_nts: a list of pointers to the Nucleotides in the Codon
        """

        self.frame = frame
        self.orf = orf
        self.nts_in_codon = ptrs_to_nts

    def nt_in_pos(self, query_nt):
        """
        Finds the position of the Nucleotide in the Codon
        :param query_nt: the Nucleotide of interest
        :return: the position of the Nucleotide in the Codon
        """
        for idx, nt in enumerate(self.nts_in_codon):
            if query_nt == nt:
                return idx

    def is_nonsyn(self, pos_in_codon, to_nt):
        """
        Finds if a substitution at the specified position results in a non-synonymous mutation
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide (A, T, G, C)
        :return: True if the substitution leads to a non-synonymous mutation,
                 False if the substitution leads to a synonymous mutation
        """
        if self.orf[0] < self.orf[1]:  # Positive strand
            codon = [str(nt) for nt in self.nts_in_codon]    # Cast all Nucleotides in the Codon to strings
        else:
            codon = [nt.complement_state for nt in self.nts_in_codon]
            to_nt = COMPLEMENT_DICT[to_nt]

        mutated_codon = codon.copy()
        mutated_codon[pos_in_codon] = to_nt

        if CODON_DICT[''.join(mutated_codon)] != CODON_DICT[''.join(codon)]:
            return True
        else:
            return False
