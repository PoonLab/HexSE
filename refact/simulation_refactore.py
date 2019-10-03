# Simulate evolution in a sequence with overlapping reading frames according to the new data structure from new_sequence_info.py
import numpy as np
import random
from new_sequence_info import Sequence


class Simulate:
    """
    Simulate evolution within a sequence throughout a phylogeny
    :ivar seq: list of nucleotides. Each nucleotide stores information about their reading frames and evolutionary rates
    :ivar tree: rooted phylogenetic tree over which seq evolves
    """

    def __init__(self, sequence , phylo_tree = None):
        self.sequence  = sequence # Object of class Sequence (Double Linked list of Nucleotides)
        self.phylo_tree = phylo_tree
        self.event_tree = self.sequence.event_tree # Tree of class EventTree (all possible mutations related to sequence)


    def get_substitution(self):
        """
        Select a substitution by moving over the event_tree according to the generation of random numbers
        """
        # Select: to nucleotide
        events = self.event_tree['total_events']
        #to_mutation = self.select_value(self.event_tree['to_nt'], events ,'events_for_nt', 'stationary_frequency')
        to_mutation = self.select_weighted_values(self.event_tree['to_nt'], events, 'events_for_nt', 'stationary_frequency')

        # Select: from nucleotide
        from_dict = self.event_tree['to_nt'][to_mutation]['from_nt']
        events_for_nt = self.event_tree['to_nt'][to_mutation]['events_for_nt']
        from_mutation = self.select_weighted_values(from_dict, events_for_nt, 'number_of_events', 'kappa')
        final_mutation = self.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]
        from_nucleotide = random.choice(final_mutation['nts_in_subs'])

        return from_nucleotide, to_mutation


    def select_weighted_values(self, dictionary, number_of_total_events, key_for_local_events, key_to_weight):
        """
        Randomly selected a key from a dictionary of events with weighted values
        :param number_of_total_events: Number of total number of events on branch
        :param number_of_local_events: key to the number of events for each specific value
        :param key_to_weight: Depending on the level of the branch, this could be:
                - the key to stationary frequency
                - the key to transition/transversion rate ratio
        """

        total_events = number_of_total_events
        temp = {}
        sum_values = 0
        # Create a temp dictionary to store the weighted values
        for key, value in dictionary.items():
            if value:
                temp[key] = (value[key_for_local_events]/total_events)*value[key_to_weight]
                sum_values += temp[key]
            else:
                temp[key] = None

        # Randomly selected a key on the dictionary
        iter_object = iter(temp.items())
        limit = random.uniform(0, sum_values)
        s = 0
        while s < limit:
            try:
                (key, value) = next(iter_object)
                if value:
                    s += value
                else:
                    pass
            except StopIteration:
                break

        return key
