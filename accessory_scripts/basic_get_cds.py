# handle = "/home/laura/Projects/ovrf/temp/iav/iav_pb1_h5n1.fa"
# pos = (24, 2295)
# with open(handle) as alignment:
#     for line in alignment:
#         if line.startswith('>'):
#             print(line, end="")
#         else:
#             cds = line[pos[0] : pos[1]]
#             print(cds)

import argparse

def get_args(parser):
    parser.add_argument(
        'file',
        help='Path to the file containing the alignment'
    )
    parser.add_argument(
        'pos',
        help='Position of the CDS to extract as startposition and end position (i.e; 3,100)'
    )

    return parser.parse_args()

def main():
    parser = argparse.ArgumentParser(
        description='Get a CDS from an alignment'
    )

    args = get_args(parser)
    handle = args.file
    pos = args.pos.split(',')
    start, end = int(pos[0]), int(pos[1])

    with open(handle) as alignment:
        for line in alignment:
            if line.startswith('>'):
                print(line, end="")
            else:
                cds = line[start : end]
                print(cds)


if __name__ == '__main__':
    main()
