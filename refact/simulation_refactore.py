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
        percentage = random.uniform(0, 1)
        total = event_tree['total_events']
        limit = total * percentage
        sum = 0

        # Draw event
        my_tree = event_tree['to_nt']
        iter_object = iter(my_tree.items())
        result = None
        while sum < limit:
            try:
                (key, to_nt) = next(iter_object)
                result = to_nt
                number = to_nt['events_for_nt']
                sum += number
                print("limit: {}, sum: {} ".format(limit, sum))
            except StopIteration:
                break

        print(result)


