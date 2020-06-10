# Get a specific reading frame from an alignment (simulation output)
import argparse
import os

handle = "/home/laura/Projects/ovrf/temp/hiv_0.05.fa"
pos = (200,500)

def get_args(parser):
    parser.add_argument(
        'file',
        help='Path to the file containing the alignment'
    )
    parser.add_argument(
        'start', type=int,
        help='Start position of the CDS'
    )
    parser.add_argument(
        'end', type=int,
        help='Final position of the CDS'
    )

    return parser.parse_args()

def get_cds(alignment_file, start, end):
    file_dir = os.path.dirname(alignment_file)
    file_name = os.path.basename(alignment_file)
    out_file = open("{}/{}_{}_{}".format(file_dir, start, end, file_name),"w")
    with open(alignment_file) as alignment:
        for line in alignment:
            if line.startswith('>'):
                out_file.write(line)
            else:
                cds = line[start : end-3]  # (-3) Remove the stop codon at the end of the reading frame
                out_file.write(cds)
                out_file.write('\n')

def main():
    parser = argparse.ArgumentParser(
        description='Get an ORF from an alignment'
    )
    args = get_args(parser)
    alignment_file = args.file
    start = args.start
    end = args.end
    get_cds(alignment_file, start, end)

if __name__ == '__main__':
    main()
