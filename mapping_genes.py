from Bio import SeqIO

# Adapted from mgalactus https://www.biostars.org/p/142972/#142977

gbk_filename = "sequence.gb"
map = {}
mapped_array = []
min_position = 99999999999999
max_position = -1

# Identify start and end position for each gene in the genbank file
for seq_record in SeqIO.parse(gbk_filename, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS":
            # print('Start: %d, Stop: %d'%(int(seq_feature.location.start),
            #                              int(seq_feature.location.end)))
            start = int(seq_feature.location.start)
            end = int(seq_feature.location.end)

            if start < min_position:
                min_position = start

            if end > max_position:
                max_position = end

            for i in range(start, end):
                if map.get(i):
                    map[i] += 1
                else:
                    map[i] = 1

# Print the overlapping ones from a map
# for key, value in map.iteritems():
#     if value > 1:
#         print key, value

# Count number of genes that map with each position
mapped_array = [0] * (max_position - min_position + 1)

for key, value in map.iteritems():
    mapped_array[key-min_position] = value

# print mapped_array

# Print the overlapping from a mapped_array
print len(mapped_array)
for i, counts in enumerate(mapped_array):
     if counts == 2:
        print i+min_position, counts
