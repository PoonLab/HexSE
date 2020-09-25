handle = "/home/laura/Projects/ovrf/temp/iav/iav_pb1_h5n1.fa"
pos = (24, 2295)
with open(handle) as alignment:
    for line in alignment:
        if line.startswith('>'):
            print(line, end="")
        else:
            cds = line[pos[0] : pos[1]]
            print(cds)
