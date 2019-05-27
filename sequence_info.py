# Store information of the sequence

def sort_orfs(orfs):
    """
    Store orfs in position according to plus zero orf (first of the list). They will be classified as (+0, +1, +2, -0, -1, -2)
    :param orfs: list of orfs as tuples for <seq> (ex. [(5,16),(11,0)])
    :return: list of ORFs classified according to their shift regarding to the plus zero one (+0, +1, +2, -0, -1, -2)
    """

    plus_cero_orf = orfs[0]
    orf_position = [1] * 6
    orf_position[0] = plus_cero_orf


    for orf in orfs[1:]:
        if type(orf) == tuple:
            if orf[0] < orf[1]: # positive strand
                difference = (orf[0] - plus_cero_orf[0]) % 3
                if difference == 1: # plus one
                    orf_position[1] = orf
                elif difference == 2: # plus two
                    orf_position[2] = orf
            elif orf[0] > orf[1]: # negative strand
                difference = (plus_cero_orf[1] - orf[0]) % 3
                if difference == 0:
                    orf_position[3] = orf # minus zero
                elif difference == 1:
                    orf_position[4] = orf # minus one
                elif difference == 2:
                    orf_position[5] = orf # minus two

    return orf_position

class Sequence(list):
    """
    List contaning information about each nucleotide in seq
    """

    def __init__(self, original_seq, sorted_orfs):
        super(Sequence, self).__init__()
        self.original_seq = original_seq
        for position in range(len(original_seq)):
            self.append(Nucleotide(original_seq[position], position, sorted_orfs))


class Nucleotide(str):
    """
    Nucleotide information (letter, position in seq and orfs to which it belongs to)
    """

    # method '__new__': https://stackoverflow.com/a/30045261/638425
    def __new__(self, letter, *args, **kwargs):
        return super(Nucleotide, self).__new__(self, letter)

    def __init__(self, letter, position, sorted_orfs):
        """
        :param letter: Nucleotide A, C, G or T
        :param position: positon in <seq>
        :param sorted_orfs: list of orf as tuples defined by user
        """
        self.position = position
        orf_presence = [False,] * 6
        for i in range(len(sorted_orfs)):
            orf = sorted_orfs[i]
            if type(orf) == tuple:
                if orf[0] < orf[1]: # positive strand
                    if self.position in range(orf[0],orf[1]):
                        orf_presence[i] = True
                elif orf[0] > orf[1]: # negativa strand
                    if self.position in range(orf[0],orf[1], -1):
                        orf_presence[i] = True
        self.orf_presence = orf_presence
