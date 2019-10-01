# Create a tree with application of mutation rates according to the landscape of possibilities

TRANSITIONS_DICT = {'A':'G', 'G':'A', 'T':'C', 'C':'T'}

class EventTree:

    def __init__(self, nt_frequencies):
        """

        :param nt_frequencies: <Dict> Frequency of nucleotides in a given sequence, with nucleotide as keys
        :param kappa: transition/transversion rate ratio
        :param omega_values: <Dict> Numeric values drawn from gamma distribution. Applied in case that a mutation is non_synonymous

        """
        self.nt_frequencies = nt_frequencies

    def get_nt_frequencies(self):
        return self.nt_frequencies

    def check_transversion(self, from_nt, to_nt):
        trv_dict = {'is_trv':True}
        if TRANSITIONS_DICT[from_nt] == to_nt:
            trv_dict['is_trv'] = False
        return trv_dict

    def create_event_tree(self):
        event_tree = {'to_nt': {'A':{},
                                'T':{},
                                'C':{},
                                'G':{}}
                        }

        for to_nt in self.nt_frequencies.keys():
            if to_nt in event_tree['to_nt'].keys():
                # Add stationary frequencies to every nucleotide in the event tree
                event_tree['to_nt'][to_nt]['stationary_frequency'] = self.nt_frequencies[to_nt]
                # Update nucleotides with possible mutations
                event_tree['to_nt'][to_nt].update([('from_nt',{'A':{},'T':{},'C':{},'G':{}})])
                # For possible mutations, check if they are a transition or a transversion
                for from_nt in event_tree['to_nt'][to_nt]['from_nt'].keys():
                    # Nucleotide cannot change to itself
                    if from_nt == to_nt:
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = None
                    else:
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt] = self.check_transversion(from_nt, to_nt)
                        # create key to a dictionary that will store information about nucleotides affected by nonsyn mutations
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_nonsyn', {})])
                        event_tree['to_nt'][to_nt]['from_nt'][from_nt].update([('is_syn', [])])

        return event_tree

# s_frequencies = {'A': 0.24, 'C': 0.24, 'T': 0.24, 'G': 0.29}
# events = EventTree(s_frequencies)
# dict = events.create_event_tree()
# print(dict['to_nt']['T']['from_nt']['C'])
# print(dict)
