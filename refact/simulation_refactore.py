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
        to_mutation = self.select_value(self.event_tree['to_nt'], events ,'events_for_nt')
        # Select: from nucleotide
        events = to_mutation[0]['events_for_nt']
        from_mutation = self.select_value(to_mutation[0]['from_nt'], events, 'number_of_events')

        return from_mutation


    def select_value(self, dictionary, number_of_events, key_local):
        """
        :param number_of_events: Number of total number of events on branch
        :param key_local: key to the number of events for each specific value
        """
        percentage = random.uniform(0, 1)
        limit = number_of_events * percentage
        iter_object = iter(dictionary.items())
        result = None
        sum = 0
        while sum < limit:
            try:
                (key, value) = next(iter_object)
                out_key = key
                result = value
                if value:
                    number = value[key_local]
                    sum += number
                else:
                    pass
            except StopIteration:
                break
        return result, out_key
