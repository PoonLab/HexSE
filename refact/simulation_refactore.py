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
        to_mutation = self.select_weighted_values(self.event_tree['to_nt'], events, 'events_for_nt', 'stationary_frequency')

        # Select: possible from nucleotides
        from_dict = self.event_tree['to_nt'][to_mutation]['from_nt']
        events_for_nt = self.event_tree['to_nt'][to_mutation]['events_for_nt']
        from_mutation = self.select_weighted_values(from_dict, events_for_nt, 'number_of_events', 'kappa')
        final_mutation = self.event_tree['to_nt'][to_mutation]['from_nt'][from_mutation]

        # List of nucleotides that are candidates to mutate
        candidate_nts = final_mutation['nts_in_subs']
        rates_list = [nt.get_mutation_rate() for nt in candidate_nts]
        nt_dict = dict(zip(candidate_nts, rates_list))

        # Select weighted nucleotide
        from_nucleotide = self.weighted_random_choice(nt_dict, sum(rates_list))

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
                # weight values
                temp[key] = (value[key_for_local_events]/total_events)*value[key_to_weight]
                sum_values += temp[key]
            else:
                temp[key] = None

        # Randomly selected a key on the dictionary
        return self.weighted_random_choice(temp, sum_values)

    @staticmethod
    def weighted_random_choice(dictionary, sum_values):
        """
        Randomly select a key on a dictionary where values correspond to the weight for the key
        :param dictionary: Dictionary to select the value from
        :param sum_values: sum all values on dict to establish the limit for the mutation
        :return: random key from dict
        """

        iter_object = iter(dictionary.items())
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

    def sum_rates(self):
        """
        Calculate the total mutation rate of sequence
        """
        total_rate = sum([nt.get_mutation_rate() for nt in iter(self.sequence.nt_sequence)])
        return total_rate


    def draw_waiting_time(self):
        """
        Draw a time at which mutation occurs according to mutation rates.
        :return:
        """
        instant_rate = self.sum_rates()
        time = np.random.exponential(scale=instant_rate)
        return time

    def mutate_on_branch(self, branch_length):
        """
        Simulate molecular evolution in sequence given a branch length
        """
        times_sum = 0

        while True:
            random_time = self.draw_waiting_time()
            times_sum += random_time
            if times_sum > branch_length:
                break

            # Draw a mutation
            mutation = self.get_substitution()
            # Replace state of nucleotide
            mutation[0].set_state(mutation[1])

            # TODO: update parameters of mutated and adjacent nucleotides according to the new state