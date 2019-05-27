# Store information of the sequence

def sort_orfs(orfs):
    """
    Store orfs in position according to plus zero orf (first of the list). They will be classified as (+0, +1, +2, -0, -1, -2)
    :param orfs: list of orfs for <seq> (ex. [(5,16),(11,0)])
    :return: list of ORFs classified according to their shift regarding to the plus zero one (+0, +1, +2, -0, -1, -2)
    """

    plus_cero_orf = orfs[0]
    orf_position = [1] * 6
    orf_position[0] = plus_cero_orf

    for orf in orfs[1:]:
        if orf[0] < orf[1]: # positive strand
            difference = (orf[0] - plus_cero_orf[0]) % 3
            print(difference)
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

class Nucleotide(str):

    def __init__(self, nucleotide, orfs = []):

        super(Nucleotide, self).__init__()
        self.orfs = orfs
