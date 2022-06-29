# From GenBank file create a YML to run on HexSE

import argparse
from Bio import SeqIO
import yaml
from random import randint

print("This is a script")

handle = "HIV_NC_001802.gb"

def get_locations(handle):
    
    open(handle)
    locations = []
    for record in SeqIO.parse(handle, format="genbank"):
        cds = [feat for feat in record.features if feat.type=="CDS"]
        for cd in cds:
            for loc in cd.location.parts:
                start, end = loc.start, loc.end
                locations.append((int(start),int(end)))

    return locations


orf_info = {}
for location in get_locations(handle):
    loc_string = ','.join([str(value) for value in location])
    omega_shape = 1+ (randint(0,10)/10)
    orf_info[loc_string] = {
                            'omega_classes': randint(2,6),
                            'omega_shape': omega_shape,
                            'omega_dist': 'gamma'
                        }

info_for_yaml = {
                    'global_rate': 0.5,
                    'kappa': 0.3,
                    'pi': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T':0.25},
                    'mu':{'classes': 2, 'shape': 1.0, 'dist': 'lognorm'},
                    'circular': 'false',
                    'orfs': orf_info
                }


with open('NC_001802_HIV.yaml', 'w') as file:
    documents = yaml.dump(info_for_yaml, file)