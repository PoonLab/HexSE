# Store sequence information

import random
import copy
import sys

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

    def __init__(self, str_sequence, orfs, kappa, mu, pi, omega_values, cat_values, circular=False):
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
        self.omega_values = omega_values  # Numeric values for omega (drawn from a gamma distribution by default)
        self.__codons = []      # Store references to all codons
        self.nt_sequence = []   # List of Nucleotide objects
        self.is_circular = circular  # True if the genome is circular, False otherwise
        self.cat_values = cat_values # Values drawn from a gamma distribution to categorize nucleotides according to their mutation rates

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
            self.set_substitution_rates(nt)  # Get substitution rates for the nucleotide
            self.nt_in_event_tree(nt)  # Locate nucleotide in the event tree

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
    def complement(seq, rev=False):
        """
        Generates the complement of a DNA sequence
        :param seq: the input sequence
        :param <option> rev: option to find the reverse complement
        :return s: The complement (or reverse complement) of the sequence
        """
        if rev:
            s = reversed(seq.upper())
        else:
            s = seq

        result = ''
        for i in s:
            result += COMPLEMENT_DICT[i]

        return result

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
        :param nt_categories: Dictionary with classes to store nucleotides with the same mutation rates
        :return event tree: a nested dictionary containing information about the mutation event
        """
        event_tree = {'to_nt': {'A': {}, 'T': {}, 'C': {}, 'G': {}}}


        for to_nt in self.pi.keys():
            if to_nt in event_tree['to_nt'].keys():
                # Update nucleotides with possible mutations
                event_tree['to_nt'][to_nt].update([('from_nt', {'A': {'class':{}}, 'T': {'class':{}}, 'C': {'class':{}}, 'G': {'class':{}}})])

                # For possible mutations, check if they are a transition or a transversion
                for from_nt in event_tree['to_nt'][to_nt]['from_nt'].keys():
                    # Nucleotide cannot change to itself
                    if from_nt == to_nt:
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = None
                    # else:
                    #     trv_dict = {'is_trv': True, 'kappa': self.kappa}
                    #     # If the mutation is a transition, set kappa to 1
                    #     if not self.is_transv(from_nt, to_nt):
                    #         trv_dict['is_trv'] = False
                    #         trv_dict['kappa'] = 1
                    #
                    #     event_tree['to_nt'][to_nt]['from_nt'][from_nt] = trv_dict
                    #
                    #     # Create key that will store information about nucleotides affected by syn and non syn mutation in the possible classes
                    #     # Deep copies are used to avoid cross referencing the newly created dictionaries
                    #     nt_categories_copy_true = copy.deepcopy(self.cat_values)
                    #     nt_categories_copy_false = copy.deepcopy(self.cat_values)
                    #     event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_nonsyn', {'True': nt_categories_copy_true, 'False':nt_categories_copy_false})])

        return event_tree


    def set_substitution_rates(self, nt):
        """
        Calculates and sets substitution rates of a nucleotide
        :param nt: object of class Nucleotide
        :return: 1. Dictionary of substitutions rates, keyed by nt subs
                 2. Dictionary with omega keys
        """
        current_nt = nt.state
        sub_rates = {}
        my_omega_keys = {}
        my_categories_keys = {}

        for to_nt in NUCLEOTIDES:

            if to_nt == current_nt:
                sub_rates[to_nt] = None
                my_omega_keys[to_nt] = None

            else:
                # Apply global substitution rate and stationary nucleotide frequency
                sub_rates[to_nt] = self.mu * self.pi[current_nt]
                if self.is_transv(current_nt, to_nt):
                    sub_rates[to_nt] *= self.kappa

                chosen_omegas = []
                # If nucleotide belongs to a codon
                if nt.codons:

                    # If mutation does not introduce a STOP and nucleotide is not part of a START codon
                    if not self.is_start_stop_codon(nt, to_nt):
                        for codon in nt.codons:
                            pos_in_codon = codon.nt_in_pos(nt)

                            # Apply omega when mutation is non-synonymous
                            if codon.is_nonsyn(pos_in_codon, to_nt):
                                selected_omega = random.choice(list(self.omega_values))  # Radomly select a key in the omega values dictionary
                                sub_rates[to_nt] *= self.omega_values[selected_omega]
                                chosen_omegas.append(selected_omega)

                    else:  # Mutation introduce a STOP or nt is part of a START
                        sub_rates[to_nt] *= 0

                    # Randomly select one of the mu Values
                    selected_cat = random.choice(list(self.cat_values))
                    sub_rates[to_nt] *=  self.cat_values[selected_cat]

                    my_omega_keys[to_nt] = chosen_omegas  # Store omega keys used in the substitution
                                                          # If mutation is synonymous, omega_keys list will be empty
                    my_categories_keys[to_nt] = selected_cat

        # Set substitution rates and key values for the nucleotide object
        nt.set_rates(sub_rates)
        nt.set_categories(my_categories_keys)
        nt.set_omega(my_omega_keys)
        nt.get_mutation_rate()
        print(nt.omega_keys)
        print(nt.categories_keys)

    @staticmethod
    def is_start_stop_codon(nt, to_nt):
        """"
        Check if mutation is a STOP codon or nucleotide belongs to a START codon
        :return: False if mutation does not create a STOP in any of the codons the nucleotide is part of
        """
        for codon in nt.codons:

            if codon.is_stop(codon.nt_in_pos(nt), to_nt) or codon.is_start():
                return True

        return False


    def nt_in_event_tree(self, nt):
        """
        Store nucleotide in each branch of the Event Tree where it belongs
        """
        current_nt = nt.state

        for to_nt in NUCLEOTIDES:

            if to_nt != current_nt and not self.is_start_stop_codon(nt, to_nt):
                omega_class = "_".join(nt.omega_keys[to_nt])  # Create one nucleotide key with all the omegas for that substitution
                category = nt.categories_keys[to_nt]
                branch = self.event_tree['to_nt'][to_nt]['from_nt'][current_nt]['class']

                if category in branch:  # If category of the nucleotide is already a key on the event tree
                    cat_branch = branch[category]  # Enter that branch

                else:  # Category key has not been created on the Event Tree
                    branch[category] = {'syn_mutations': []}  # Create new branch
                    cat_branch = branch[category]

                # Store nucleotide according to omega keys
                if omega_class:  # substitution is nonsyn for at least one codon
                    if omega_class in cat_branch:
                        cat_branch[omega_class].append(nt)
                    else:
                        cat_branch[omega_class] = [nt]

                else:  # Mutation is syn in all Codons
                    cat_branch['syn_mutations'].append(nt)



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

    def find_codons(self, frame, orf_coords):
        """
        Gets the Codon sequence
        :param frame: the frame of the ORF
        :param orf: tuple containing the coordinates of the ORF
        :return: a list of Codon objects for the specified ORF
        """
        codons = []
        cds = []
        for coord in orf_coords:
            cds.extend(self.nt_sequence[coord[0]: coord[1]])

        # Reverse strand orf
        if frame.startswith('-'):
            cds = cds[::-1]

        # Iterate over list by threes and create Codons
        for i in range(3, len(cds) + 1, 3):
            cdn = cds[i - 3: i]
            codon = Codon(frame, orf_coords, cdn)
            codons.append(codon)

        return codons

    def check_event_tree(self):
        """
        When debugging, useful to check if nucleotides are being properly stored on the Event Tree
        """
        for key1, to_nt in self.event_tree['to_nt'].items():
            subset = to_nt['from_nt']

            for key2, from_nt in subset.items():
                if key2 != 'T' and from_nt and from_nt.get('nts_in_subs'):
                    nts_in_subs = list(from_nt['nts_in_subs'].keys())
                    if len([1 for tip in nts_in_subs if str(tip).lower() == 't0']) > 0:
                        meta2 = {'nts_in_subs': nts_in_subs}
                        print(f'>>>>>>>>>>>> meta2: from {key2}', meta2)
                        sys.exit(1)


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
        self.omega_keys = {}  # omega keys chosen when calculating rates
        self.categories_keys = {}  # category keys chosen when calculating rates
        self.mutation_rate = 0  # The total mutation rate
        self.relevant_info = {}

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
        new_nucletotide.mutation_rate = copy.deepcopy(self.mutation_rate, memodict)
        new_nucletotide.codons = []     # References to Codons will be set when the Sequence is deep-copied

        return new_nucletotide

    def set_relevant_info(self, relevant_info):
        self.relevant_info = relevant_info

    def set_state(self, new_state):
        self.state = new_state

    def get_complement_state(self):
        return self.complement_state

    def set_complement_state(self):
        self.complement_state = COMPLEMENT_DICT[self.state]

    def set_rates(self, rates):
        self.rates = rates

    def set_omega(self, omega_keys):
        self.omega_keys = omega_keys

    def set_categories(self, cat_keys):
        self.categories_keys = cat_keys

    def add_codon(self, codon):
        self.codons.append(codon)

    def get_mutation_rate(self):
        total_rate = 0
        for to_nt, value in self.rates.items():
            if value:
                total_rate += value
        self.mutation_rate = total_rate

    def get_relevant_info(self):
        """
        Create a dictionary with all relevant information related with the nucleotide
        (Useful for debugging)
        """
        info = {"state": self.state, "position": self.pos_in_seq,
         "dn": self.dN_keys, "ds":self.dS_keys,
        "rates": self.rates, "mutation rate" :self.mutation_rate, "codons": self.codons}

        return info


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

    def mutate_codon(self, pos_in_codon, to_nt):
        """
        Changes the state of the specified nucleotide in the codon
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide
        :return codon, mutated_codon: the codon and mutated codon represented as lists of strings
        """
        # Positive strand
        if self.frame.startswith('+'):
            codon = [str(nt) for nt in self.nts_in_codon]  # Cast all Nucleotides in the Codon to strings
        else:
            codon = [nt.complement_state for nt in self.nts_in_codon]
            to_nt = COMPLEMENT_DICT[to_nt]

        mutated_codon = codon.copy()
        mutated_codon[pos_in_codon] = to_nt

        return codon, mutated_codon

    def is_nonsyn(self, pos_in_codon, to_nt):
        """
        Finds if a substitution at the specified position results in a non-synonymous mutation
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide (A, T, G, C)
        :return: True if the substitution leads to a non-synonymous mutation,
                 False if the substitution leads to a synonymous mutation
        """
        codon, mutated_codon = self.mutate_codon(pos_in_codon, to_nt)
        return CODON_DICT[''.join(mutated_codon)] != CODON_DICT[''.join(codon)]

    def is_stop(self, pos_in_codon, to_nt):
        """
        Finds if a substitution at the specified position results in a stop codon
        :param pos_in_codon: the position in the Codon
        :param to_nt: the new state of the Nucleotide (A, T, G, C)
        :return: True if the substitution leads to stop codon,
                 False if the substitution doesn't lead to a stop codon
        """
        codon, mutated_codon = self.mutate_codon(pos_in_codon, to_nt)
        return CODON_DICT[''.join(mutated_codon)] == "*"

    def is_start(self):
        """
        Checks if the codon is a start codon
        :return True of the codon is a start codon, False otherwise
        """
        # Postive strand
        if self.frame.startswith('+'):
            codon = ''.join(str(nt) for nt in self.nts_in_codon)  # Cast all Nucleotides in the Codon to strings
        else:
            codon = ''.join(nt.complement_state for nt in self.nts_in_codon)

        return codon == 'ATG' and self.nts_in_codon[0].pos_in_seq == self.orf[0][0]
