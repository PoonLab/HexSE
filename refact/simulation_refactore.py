# Simulate evolution in a sequence with overlapping reading frames according to the new data structure from new_sequence_info.py

class Simulate:
    """
    Simulate evolution within a sequence throughout a phylogeny
    :ivar seq: list of nucleotides. Each nucleotide stores information about their reading frames and evolutionary rates
    :ivar tree: rooted phylogenetic tree over which seq evolves
    """

    def __init__(self, seq, tree):
        self.seq = seq
        self.tree = tree

        def sum_rates(self):
