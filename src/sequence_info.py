# Store sequence information

import copy
import random

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
    """

    def __init__(self, str_sequence, orfs, kappa, mu, pi, omegas):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.
        :param sorted_orfs: A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        :param kappa: transition/ transversion rate ratio
        :param mu: The global rate (substitutions/site/unit time)
        :param pi: Frequency of nucleotides in a given sequence, with nucleotide as keys
        :param omegas: <Dict> Numeric values drawn from gamma distribution.
                       Applied in case that a mutation is non_synonymous
        """
        self.orfs = orfs  # Dictionary of of ORFs sorted by reading frame
        self.kappa = kappa  # Transition/ transversion rate ratio
        self.mu = mu  # The global rate (substitutions/site/unit time)
        self.pi = pi  # Frequency of nucleotides, with nucleotide as keys
        self.omegas = omegas  # Numeric values drawn from gamma distribution

        # Create Nucleotides
        self.nt_sequence = DoubleLinkedList()  # A List of Nucleotide objects
        for pos_in_seq, nt in enumerate(str_sequence):
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

        # Create event tree containing all possible mutations and the parameters needed to calculate the rates
        self.event_tree = self.create_event_tree()  # Nested dict containing info about all possible mutation events

        # Calculate mutation rates for each nucleotide in sequence, populate the event tree which each nucleotide
        for nt in iter(self.nt_sequence):
            sub_rates = self.get_substitution_rates(nt)
            nt.set_rates(sub_rates[0])
            nt.set_my_omegas(sub_rates[1])
            nt.get_mutation_rate()

        # Update event_tree to include a list of nucleotides on the tips
        self.event_tree = self.get_nts_on_tips()

    def __deepcopy__(self, memo):

        return self.__class__(
            str_sequence=copy.deepcopy(self.get_string_sequence(), memo),
            orfs=copy.deepcopy(self.orfs, memo),
            kappa=copy.deepcopy(self.kappa, memo),
            mu=copy.deepcopy(self.mu, memo),
            pi=copy.deepcopy(self.pi, memo),
            omegas=copy.deepcopy(self.omegas, memo))

    def get_sequence(self):
        return self.nt_sequence

    def get_event_tree(self):
        return self.event_tree

    def get_string_sequence(self):
        seq_as_string = ''
        for nt in iter(self.nt_sequence):
            seq_as_string += nt.state
        return seq_as_string

    @staticmethod
    def get_frequency_rates(seq):
        """
        Frequency of nucleotides in the DNA sequence
        :param seq: the DNA sequence
        :return frequencies: a dictionary frequencies where the key is the nucleotide and the value is the frequency
        """
        frequencies = {'A': 0, 'C': 0, 'T': 0, 'G': 0}

        for nucleotide in frequencies:
            frequencies[nucleotide] = round((float(seq.count(nucleotide)) / (len(seq))), 2)

        return frequencies

    def create_event_tree(self):
        """
        Create an event tree (nested dictionaries) that stores pi, kappa, mu,
        and information about whether a mutation is a transition or transversion.
        :return event tree: a nested dictionary containing information about the mutation event
        """
        event_tree = {'to_nt': {'A': {}, 'T': {}, 'C': {}, 'G': {}}}

        for to_nt in self.pi.keys():
            if to_nt in event_tree['to_nt'].keys():
                # Add stationary frequencies to every nucleotide in the event tree
                event_tree['to_nt'][to_nt]['stationary_frequency'] = self.pi[to_nt]
                # Update nucleotides with possible mutations
                event_tree['to_nt'][to_nt].update([('from_nt', {'A': {}, 'T': {}, 'C': {}, 'G': {}})])

                # For possible mutations, check if they are a transition or a transversion
                for from_nt in event_tree['to_nt'][to_nt]['from_nt'].keys():
                    # Nucleotide cannot change to itself
                    if from_nt == to_nt:
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = None
                    else:
                        trv_dict = {'is_trv': True, 'kappa': self.kappa}
                        # If the mutation is a transition, set kappa to 1
                        if not self.is_transv(from_nt, to_nt):
                            trv_dict['is_trv'] = False
                            trv_dict['kappa'] = 1

                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = trv_dict

                        # Create key that will store information about nucleotides affected by nonsyn mutations
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_nonsyn', {})])
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_syn', [])])

        return event_tree

    def get_nts_on_tips(self):
        """
        Look at the tips of the event tree and create keys and values for nucleotides in substitution and number of events
        :return: update_event_tree: A nested dictionary containing:
                                    - List of nucleotides for each tip
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
                    nt_subs_length = len(nt_in_substitution)
                    non_syn_subs = from_nt['is_nonsyn']
                    for key3, nts in non_syn_subs.items():
                        nt_in_substitution.extend([nt for nt in nts])

                    updated_event_tree['to_nt'][key1]['from_nt'][key2].update([('nts_in_subs', nt_in_substitution)])
                    updated_event_tree['to_nt'][key1]['from_nt'][key2].update([('number_of_events', nt_subs_length)])
                    events_for_to_nt += nt_subs_length
                    total_events += nt_subs_length

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
        current_nt = nt.state

        for to_nt in NUCLEOTIDES:
            if to_nt == current_nt:
                sub_rates[to_nt] = None
                my_omega_keys[to_nt] = None
            else:
                # Apply global substitution rate and stationary nucleotide frequency
                sub_rates[to_nt] = self.mu * self.pi[current_nt]
                if self.is_transv(current_nt, to_nt):
                    sub_rates[to_nt] *= self.kappa

                chosen_omegas = [0, 0, 0, 0]  # omegas applied given a substitution from current_nt to to_nt
                for codon in nt.codons:
                    pos_in_codon = codon.nt_in_pos(nt)
                    if codon.is_stop(pos_in_codon, to_nt): # If mutation leads to a stop codon
                        sub_rates[to_nt] *= 0
                    elif codon.is_nonsyn(pos_in_codon, to_nt):  # Apply omega when mutation is non-synonymous
                        omega_index = random.randrange(len(self.omegas))
                        sub_rates[to_nt] *= self.omegas[omega_index]
                        chosen_omegas[omega_index] += 1

                omegas_in_sub = tuple(chosen_omegas)
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
        """
        Checks if a mutation is a transition or a transversion
        :param from_nt: the current nucleotide
        :param to_nt: the new nucleotide
        :return transv: True if the mutation is a transversion,
                        False if the mutation is a transition,
                        None if the current and new nucleotides are the same
        """
        if from_nt == to_nt:
            transv = None
        else:
            transv = True
            if TRANSITIONS_DICT[from_nt] == to_nt:
                transv = False
        return transv

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
            yield my_orf[i:i + 3]
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
        my_orf = self.nt_sequence.slice_sequence(start_pos, end_pos-1)

        # Iterate over list by threes and create Codons
        for cdn in self.codon_iterator(my_orf, start_pos, end_pos-1):
            codon = Codon(frame, orf, cdn)
            codons.append(codon)

        return codons


class Nucleotide:
    """
    Stores information about the base, the open reading frames to which a Nucleotide belongs,
    and references to the previous and next base in the sequence.
    """

    def __init__(self, state, pos_in_seq, left_nt=None, right_nt=None):
        """
        :param state: Nucleotide A, C, G or T
        :param pos_in_seq: Position of the nucleotide in the sequence
        :param left_nt : reference to the adjacent nucleotide to the left (default to None)
        :param right_nt: reference to the adjacent nucleotide to the right (default to None)
        """
        self.state = state  # The nucleotide base (A, T, G, C)
        self.pos_in_seq = pos_in_seq  # The nt's position relative to the start of the sequence
        self.codons = []  # A list of codon objects the Nucleotide is part of
        self.complement_state = COMPLEMENT_DICT[self.state]  # The complement state
        self.rates = {}  # A dictionary of mutation rates
        self.my_omegas = []  # Omegas chosen when calculating rates
        self.mutation_rate = 0  # The total mutation rate
        self.left_nt = right_nt  # Reference to the nucleotide's left neighbour (prev)
        self.right_nt = left_nt  # Reference to the nucleotide's right neighbour (next)

    def __str__(self):
        return self.state

    def __repr__(self):
        return self.state.lower() + str(self.pos_in_seq)

    def __deepcopy__(self, memo):
        return self.__class__(state=copy.deepcopy(self.state, memo),
                              pos_in_seq=copy.deepcopy(self.pos_in_seq, memo))

    def set_state(self, new_state):
        self.state = new_state

    def set_pos_in_seq(self, new_pos):
        self.pos_in_seq = new_pos

    def set_left_nt(self, new_left_nt):
        self.left_nt = new_left_nt

    def set_right_nt(self, new_right_nt):
        self.right_nt = new_right_nt

    def get_complement_state(self):
        return self.complement_state

    def set_complement_state(self):
        self.complement_state = COMPLEMENT_DICT[self.state]

    def set_rates(self, rates):
        self.rates = rates

    def set_my_omegas(self, omegas):
        self.my_omegas = omegas

    def add_codon(self, codon):
        self.codons.append(codon)

    def get_mutation_rate(self):
        total_rate = 0
        for to_nt, value in self.rates.items():
            if value:
                total_rate += value
        self.mutation_rate = total_rate


class DoubleLinkedList:
    """
    Circular Double linked list linking together objects of class Nucleotide
    Default initialization with empty head node
    """

    # Constructor for circular double-linked list
    def __init__(self):
        self.head = None        # First nucleotide in the sequence
        self.tail = None        # Last nucleotide in the sequence

    # Iterator for circular double-linked list
    def __iter__(self):
        curr_nt = self.head
        while curr_nt is not None:
            yield curr_nt
            curr_nt = curr_nt.right_nt

    # Iterator for circular double-linked list
    def __next__(self):
        curr_nt = self.head
        if curr_nt is None:
            raise StopIteration
        else:
            curr_nt = curr_nt.right_nt
            return curr_nt

    def __deepcopy__(self):
        return self.__class__()

    # Insert for circular double linked list
    def insert_nt(self, state, position):
        """
        Insert objects of class Nucleotide to the end of the DoubleLinkedList
        :param state: Nucleotide state in sequence
        :param position: Position of nt in sequence
        """
        new_nt = Nucleotide(state, position)

        # If the circular double linked list is empty, assign the first nucleotide as the head
        if self.head is None:
            self.head = new_nt
            new_nt.left_nt = None
            new_nt.right_nt = None
            self.tail = new_nt

        else:
            # Traverse the list to find the last nucleotide
            curr_nt = self.head
            while curr_nt.right_nt is not None:
                curr_nt = curr_nt.right_nt
            # Update the references
            curr_nt.right_nt = new_nt
            new_nt.left_nt = curr_nt
            self.tail = new_nt

    def nucleotide_at_pos(self, position):
        """
        Traverse sequence to find the Nucleotide object at a specific position
        :param position: the position of the Nucleotide
        :return: the Nucleotide object in the specified position
        """
        current_nt = self.head
        while current_nt is not None:
            if position == current_nt.pos_in_seq:
                return current_nt
            current_nt = current_nt.right_nt
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
        else:  # Negative strand
            curr_nt = self.nucleotide_at_pos(end_pos)
            end_pos = start_pos

        while curr_nt is not None and curr_nt.pos_in_seq <= end_pos:
            sub_seq.append(curr_nt)
            curr_nt = curr_nt.right_nt

        return sub_seq

    # Constructor for double linked list
    # def __init__(self):
    #     self.head = None            # The starting nucleotide
    #     self.current_nt = None      # Pointer to the current nucleotide for insertion
    #     self.next_iter_nt = None    # The current state of the iteration
    #     self.size = 0

    # Iterator for double linked list
    # def __iter__(self):
    #     self.next_iter_nt = self.head
    #     return self

    # Iterator for double linked list
    # def __next__(self):
    #     # Note: Not thread safe
    #     if self.next_iter_nt is not None:
    #         nt = self.next_iter_nt
    #         self.next_iter_nt = nt.right_nt
    #         return nt
    #     else:
    #         raise StopIteration

    # Insert for double linked list
    # def insert_nt(self, state, position):
    #     """
    #     Insert objects of class Nucleotide to the end of the DoubleLinkedList
    #     :param state: Nucleotide state in sequence
    #     :param position: Position of nt in sequence
    #     """
    #     new_nt = Nucleotide(state, position)  # create new Nucleotide object
    #     self.size += 1
    #
    #     # Assign the first nucleotide as head
    #     if self.head is None:
    #         self.head = new_nt
    #         self.current_nt = new_nt
    #     else:
    #         new_nt.set_left_nt(self.current_nt)  # For the new nucleotide, create a left pointer towards the current one
    #         self.current_nt.set_right_nt(new_nt)  # Create the double link between current and new
    #         self.current_nt = new_nt


class Codon:
    """
    Stores information about the frameshift, ORF, and pointers to 3 Nucleotide objects
    """

    def __init__(self, frame, orf, nts_in_codon):
        """
        Create a Codon
        :param frame: the reading frame (+0, +1, +2, -0, -1, -2)
        :param orf: a tuple containing the reading frame and the coordinates of the orf
        :param nts_in_codon: a list of pointers to the Nucleotides in the Codon
        """
        self.frame = frame  # The reading frame
        self.orf = orf  # Tuple containing the reading frame and the coordinates
        self.nts_in_codon = nts_in_codon  # List of Nucleotides in the Codon

    def __deepcopy__(self, memo):
        return self.__class__(
            frame=copy.deepcopy(self.frame, memo),
            orf=copy.deepcopy(self.orf, memo),
            nts_in_codon=copy.deepcopy(self.nts_in_codon, memo))

    def nt_in_pos(self, query_nt):
        """
        Finds the position of the Nucleotide in the Codon
        :param query_nt: the Nucleotide of interest
        :return: the position of the Nucleotide in the Codon
        """
        for idx, nt in enumerate(self.nts_in_codon):
            if query_nt is nt:
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
            codon = [str(nt) for nt in self.nts_in_codon]  # Cast all Nucleotides in the Codon to strings
        else:
            codon = [nt.complement_state for nt in self.nts_in_codon]
            to_nt = COMPLEMENT_DICT[to_nt]

        mutated_codon = codon.copy()
        mutated_codon[pos_in_codon] = to_nt

        if CODON_DICT[''.join(mutated_codon)] != CODON_DICT[''.join(codon)]:
            return True
        else:
            return False

    def is_stop(self, pos_in_codon, to_nt):
            """
            Finds if a substitution at the specified position results in a stop codon
            :param pos_in_codon: the position in the Codon
            :param to_nt: the new state of the Nucleotide (A, T, G, C)
            :return: True if the substitution leads to stop codon,
                     False if the substitution doesn't lead to a stop codon
            """
            if self.orf[0] < self.orf[1]:  # Positive strand
                codon = [str(nt) for nt in self.nts_in_codon]  # Cast all Nucleotides in the Codon to strings
            else:
                codon = [nt.complement_state for nt in self.nts_in_codon]
                to_nt = COMPLEMENT_DICT[to_nt]

            mutated_codon = codon.copy()
            mutated_codon[pos_in_codon] = to_nt

            if CODON_DICT[''.join(mutated_codon)] == "*":
                return True
            else:
                #print(''.join(mutated_codon))
                return False
