from Bio import SeqIO
handle = "HIV.gb"
#record_seq = SeqIO.parse(handle, format = "genbank")

# def retrieve_coord(handle):
#     orfs = []
#     # Loop trough records.
#     for rec in SeqIO.parse(handle, format = "genbank"):
#         full_seq = rec.seq  #  TO DO: deal with multipartite viruses?
#         cds = [feat for feat in rec.features if feat.type=="CDS"]
#         for cd in cds:
#             coord = ([(int(loc.start), int(loc.end)) for loc in cd.location.parts])  # Extract coordinates
#             print("coordinates", coord)
#             orfs.extend(coord)
#             print(cd.extract(seq))
#     return orfs, full_seq

# handle = handle.lower()
# if handle.endswith(".gb") or handle.endswith("genbank"):
#     orfs, seq = retrieve_coord(handle)
#     # print(orfs)
#     print((orfs))
#     #print(seq[265:13468])
#     # print(retrieve_coord(handle))
# elif handle.endswith(".fasta") or handle.endswith(".fa"):
#     print("This is a fasta file")
# else:
#     print("Unrecognized format")


# def retrieve_

def parse_genbank_orfs(seq_path): #TODO: change name of input
    """
    Extract ORFs from the GenBank file
    """
    orf_locations = {'+': [], '-': []}  # ORF locations sorted by strand

    # Loop through records
    for rec in SeqIO.parse(seq_path, format="genbank"):
        # Read ORFs from GenBank file
        cds = [feat for feat in rec.features if feat.type == "CDS"]
        # Record the first occurrence of the ORFs
        for cd in cds:
            orf = {'coords': []}
            strand = ''

            for loc in cd.location.parts:
                # Get the strand
                if loc.strand > 0:
                    strand = '+'
                else:
                    strand = '-'
                orf['coords'].append((int(loc.start), int(loc.end)))

            orf_locations[strand].append(orf)

    return orf_locations

print("I am running")
print(parse_genbank_orfs(handle))