from Bio import SeqIO
handle = "../temp/hiv.gb"
#record_seq = SeqIO.parse(handle, format = "genbank")

def retrieve_coord(handle):
    orfs = []
    # Loop trough records.
    for rec in SeqIO.parse(handle, format = "genbank"):
        full_seq = rec.seq  #  TO DO: deal with multipartite viruses?
        cds = [feat for feat in rec.features if feat.type=="CDS"]
        for cd in cds:
            coord = ([(int(loc.start), int(loc.end)) for loc in cd.location.parts])  # Extract coordinates
            print("coordinates", coord)
            orfs.extend(coord)
            print(cd.extract(seq))
    return orfs, full_seq

handle = handle.lower()
if handle.endswith(".gb") or handle.endswith("genbank"):
    orfs, seq = retrieve_coord(handle)
    # print(orfs)
    print((orfs))
    #print(seq[265:13468])
    # print(retrieve_coord(handle))
elif handle.endswith(".fasta") or handle.endswith(".fa"):
    print("This is a fasta file")
else:
    print("Unrecognized format")


def retrieve_
