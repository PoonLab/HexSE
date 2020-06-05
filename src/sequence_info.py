# Store sequence information

import random
import copy

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

    def __init__(self, str_sequence, orfs, kappa, mu, pi, dN_values, dS_values, circular=False):
        """
        Creates a list of nucleotides, locates open reading frames, and creates a list of codons.
        :param orfs: A dictionary of ORFs, sorted by reading frame where:
                        - the keys are the reading frames (+0, +1, +2, -0, -1, -2)
                        - the values are a list of tuples containing the start and end positions of the ORFs.
                        - Ex: {'+0': [(0, 8), (3, 15)], '+1': [(1, 9)], '+2': [], '-0': [], '-1': [], '-2': []}
        :param kappa: transition/ transversion rate ratio
        :param mu: The global rate (substitutions/site/unit time)
        :param pi: Frequency of nucleotides in a given sequence, with nucleotide as keys
        :param omegas: <Dict> Numeric values drawn from gamma distribution.
                       Applied in case that a mutation is non_synonymous
        :param circular: True if the genome is circular, false if the genome is linear (default: linear)
        """
        self.orfs = orfs  # Dictionary of of ORFs sorted by reading frame
        self.kappa = kappa  # Transition/ transversion rate ratio
        self.mu = mu  # The global rate (substitutions/site/unit time)
        self.pi = pi  # Frequency of nucleotides, with nucleotide as keys
        self.dN_values = dN_values  # Numeric values for dN (drawn from a gamma distribution by default)
        self.dS_values = dS_values  # Numeric values for dS (drawn from a gamma distribution by default)
        self.__codons = []      # Store references to all codons
        self.nt_sequence = []   # List of Nucleotide objects
        self.is_circular = circular  # True if the genome is circular, False otherwise

        # Create Nucleotides
        for pos_in_seq, nt in enumerate(str_sequence):
            self.nt_sequence.append(Nucleotide(nt, pos_in_seq))

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
                        self.__codons.append(codon)

        # Create event tree containing all possible mutations and the parameters needed to calculate the rates
        self.event_tree = self.create_event_tree()  # Nested dict containing info about all possible mutation events

        # Calculate mutation rates for each nucleotide in sequence, populate the event tree which each nucleotide
        for nt in self.nt_sequence:
            sub_rates = self.get_substitution_rates(nt)
            nt.set_rates(sub_rates[0])
            nt.set_dN(sub_rates[1])
            nt.set_dS(sub_rates[2])
            nt.get_mutation_rate()

        # Update event_tree to include a list of nucleotides on the tips
        self.event_tree = self.get_nts_on_tips()

    def __deepcopy__(self, memodict):
        """
        Creates a deepcopy of Sequence and sets the reference(s) for a Nucleotide's Codon(s)
        """

        # Creates a new Sequence
        cls = self.__class__
        new_sequence = cls.__new__(cls)
        memodict[id(self)] = new_sequence        # Avoid duplicate copying

        # Set attributes of new Sequence to the same as the original object
        for k, v in self.__dict__.items():
            setattr(new_sequence, k, copy.deepcopy(v, memodict))

        # Set references to Codons
        for codon in new_sequence.__codons:
            for i, nt in enumerate(codon.nts_in_codon):
                nt.codons.append(codon)

        return new_sequence

    def __str__(self):
        """
        Represents the Sequence as a string by casting nt_sequence to a string.
        """
        return ''.join(str(nt) for nt in self.nt_sequence)

    def get_sequence(self):
        return self.nt_sequence

    def get_right_nt(self, pos_in_seq):
        """
        Returns the next Nucleotide in the sequence
        :param pos_in_seq: the position of the Nucleotide in the sequence
        """
        if self.is_circular:
            if pos_in_seq == len(self.nt_sequence) - 1:
                return self.nt_sequence[0]
        else:
            return self.nt_sequence[pos_in_seq + 1]

    def get_left_nt(self, pos_in_seq):
        """
        Returns the previous Nucleotide in the sequence
        :param pos_in_seq: the position of the Nucleotide in the sequence
        """
        if pos_in_seq == 0:
            if self.is_circular:
                return self.nt_sequence[- 1]
            else:
                return None
        else:
            return self.nt_sequence[pos_in_seq - 1]

    def get_event_tree(self):
        return self.event_tree

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
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_nonsyn', {'dN': {}, 'dS': {}})])
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
                    non_syn_subs = from_nt['is_nonsyn']['dN']
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
        my_dN_keys = {}
        my_dS_keys = {}
        current_nt = nt.state

        for to_nt in NUCLEOTIDES:
            if to_nt == current_nt:
                sub_rates[to_nt] = None
                my_dN_keys[to_nt] = None
                my_dS_keys[to_nt] = None
            else:
                # Apply global substitution rate and stationary nucleotide frequency
                sub_rates[to_nt] = self.mu * self.pi[current_nt]
                if self.is_transv(current_nt, to_nt):
                    sub_rates[to_nt] *= self.kappa

                # dN and dS applied given a substitution from current_nt to to_nt
                chosen_dN = [0 for _ in range(len(self.dN_values))]
                chosen_dS = [0 for _ in range(len(self.dS_values))]

                for codon in nt.codons:
                    pos_in_codon = codon.nt_in_pos(nt)
                    if codon.is_stop(pos_in_codon, to_nt):  # If mutation leads to a stop codon
                        sub_rates[to_nt] *= 0
                        nt.is_stop = True
                    elif codon.is_nonsyn(pos_in_codon, to_nt):  # Apply omega when mutation is non-synonymous
                        dN_index = random.randrange(len(self.dN_values))
                        dS_index = random.randrange(len(self.dS_values))
                        sub_rates[to_nt] *= (self.dN_values[dN_index] / self.dS_values[dS_index])
                        chosen_dN[dN_index] += 1
                        chosen_dS[dS_index] += 1

                dN_in_sub = tuple(chosen_dN)
                dS_in_sub = tuple(chosen_dS)
                my_dN_keys[to_nt] = dN_in_sub
                my_dS_keys[to_nt] = dS_in_sub

                current_dS = self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']['dS']
                current_dN = self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']['dN']

                # Populate event tree using dN and dS keys
                if any(dN_in_sub) and any(dS_in_sub):
                    # Create dN key if needed, associate it with the current nucleotide
                    if dN_in_sub not in current_dN:
                        self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']['dN'][dN_in_sub] = [nt]
                    else:
                        self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']['dN'][dN_in_sub].append(nt)

                    # Create dS key if needed, associate it with the current nucleotide
                    if dS_in_sub not in current_dS:
                        self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']['dS'][dS_in_sub] = [nt]
                    else:
                        self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_nonsyn']['dS'][dS_in_sub].append(nt)

                elif 0 not in sub_rates.values():  # If does not introduce an STOP codon
                    self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['is_syn'].append(nt)  # Mutation is synonymous

        return sub_rates, my_dN_keys, my_dS_keys

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

        # Reverse stand ORF
        if start_pos > end_pos:
            my_orf = self.nt_sequence[end_pos: start_pos]
        else:
            my_orf = self.nt_sequence[start_pos: end_pos]

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

    def __init__(self, state, pos_in_seq):
        """
        :param state: Nucleotide A, C, G or T
        :param pos_in_seq: Position of the nucleotide in the sequence
        """
        self.state = state  # The nucleotide base (A, T, G, C)
        self.pos_in_seq = pos_in_seq  # The nt's position relative to the start of the sequence
        self.codons = []  # A list of codon objects the Nucleotide is part of
        self.complement_state = COMPLEMENT_DICT[self.state]  # The complement state
        self.rates = {}  # A dictionary of mutation rates
        self.dN_values = []  # dN chosen when calculating rates
        self.dS_values = []  # dS chosen when calculating rates
        self.mutation_rate = 0  # The total mutation rate

    def __str__(self):
        return self.state

    def __repr__(self):
        return self.state.lower() + str(self.pos_in_seq)

    def __deepcopy__(self, memodict):
        """
        Creates a deepcopy of a Nucleotide.
        Note: A Nucleotide's reference(s) to its Codon(s) will be set in Sequence's deepcopy
        """

        # Creates a new Nucleotide
        cls = self.__class__
        new_nucletotide = cls.__new__(cls)
        memodict[id(self)] = new_nucletotide        # Avoid duplicate copying

        # Copy all attributes except the codons
        new_nucletotide.state = copy.deepcopy(self.state, memodict)
        new_nucletotide.pos_in_seq = copy.deepcopy(self.pos_in_seq, memodict)
        new_nucletotide.complement_state = copy.deepcopy(self.complement_state, memodict)
        new_nucletotide.rates = copy.deepcopy(self.rates, memodict)
        new_nucletotide.dN_values = copy.deepcopy(self.dN_values, memodict)
        new_nucletotide.dS_values = copy.deepcopy(self.dS_values, memodict)
        new_nucletotide.mutation_rate = copy.deepcopy(self.mutation_rate, memodict)
        new_nucletotide.codons = []     # References to Codons will be set when the Sequence is deep-copied

        return new_nucletotide

    def set_state(self, new_state):
        self.state = new_state

    def get_complement_state(self):
        return self.complement_state

    def set_complement_state(self):
        self.complement_state = COMPLEMENT_DICT[self.state]

    def set_rates(self, rates):
        self.rates = rates

    def set_dN(self, dN_values):
        self.dN_values = dN_values

    def set_dS(self, dS_values):
        self.dS_values =dS_values

    def add_codon(self, codon):
        self.codons.append(codon)

    def get_mutation_rate(self):
        total_rate = 0
        for to_nt, value in self.rates.items():
            if value:
                total_rate += value
        self.mutation_rate = total_rate


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

    def __repr__(self):
        return ''.join(str(nt) for nt in self.nts_in_codon)

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
            return False
