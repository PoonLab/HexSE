# Script to extract open reading frames from a genbank file (according to CDS)
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import re

file = "NC_003977.2_HBV.gb"
#rec = SeqIO.parse(open(file,"r"), "genbank")

def orfs_from_gb(file):
    for gb_record in SeqIO.parse(open(file,"r"), "genbank") :
        # now do something with the record
        location_out = []
        for feature in gb_record.features:
            if feature.type == "CDS":
                loc = feature.location
                start = feature.location.start.position
                end = feature.location.end.position
                out = (start,end)
                location_out.append(out)
    return location_out


def sq_from_gb(file):
    for gb_record in SeqIO.parse(open(file, "r"), "genbank"):
        seq = gb_record.seq
    return seq