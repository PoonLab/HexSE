# From GenBank file create a YML to run on HexSE

import argparse
from Bio import SeqIO
import yaml
from random import randint


def get_args(parser):

    parser.add_argument(
        'gb_file',
        help = 'Path to the genbank file'
    )
    parser.add_argument(
        'yaml_file',
        help = 'Path to the yaml configuration file'
    )
    parser.add_argument(
        '--k', type=float,
        help = 'kappa: transition/transversion rate ratio',
        default = 0.3

    )
    parser.add_argument(
        '--p', type=dict,
        help = 'pi: nucleotide frequency',
        default = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T':0.25}

    )

    parser.add_argument(
        '--m',
        help = 'mu: mutation rate. Interger if only one, dictionary if values to be drawn from distribution',
        default = 1,
        # {'classes': 2, 'shape': 1.0, 'dist': 'lognorm'},

    )

    return parser.parse_args()

def get_locations(handle):
    
    open(handle)
    locations = []
    for record in SeqIO.parse(handle, format="genbank"):
        cds = [feat for feat in record.features if feat.type=="CDS"]
        for cd in cds:
            temp = []
            for loc in cd.location.parts:
                start, end = loc.start, loc.end
                temp.append((int(start),int(end)))
            locations.append(temp)

    return locations

def generate_param_values(class_range=(2,6), shape_range=(0.1,1.5), scale_range=(0.1,1.5)):
    """
    Randomly selects an interger for the number of categories and floats for shape and scale
    :param class_range: tuple, minimum and maximum values for number of classes
    :param class_shape: tuple, minimum and maximum values for shape parameter
    :param class_scale: tuple, minimum and maximum values for rate parameter
    """
    number_clases = randint(class_range[0], class_range[1])
    flt=(randint(0,10)/10)
    # Transform to get range:((max - min)float + min)
    shape = round(((shape_range[1] - shape_range[0]) * flt + shape_range[0]),1)
    scale = round(((scale_range[1] - scale_range[0]) * flt + scale_range[0]),1)

    return number_clases, shape, scale

def create_orf_dict(orf_loc, dn_values, ds_values, dist='gamma'):
    """
    Store orf parameter values in a dictionary with keys as required by HexSE

    :param orf_loc: str, orf coordinates in the genome. 
    :param dn_values: list, number of classes, shape value, and scale value to inform dn distribution
    :param ds_values: list, number of classes, shape value, and scale value to inform ds distribution
    :param: str, distribution from which dn and ds values will be drawn
    """
    
    orf_dict = {}
    orf_dict[orf_loc] = {   # dN (non-synonymous) parameters
                            'dn_class': dn_values[0],
                            'dn_dist': dist,
                            'dn_shape': dn_values[1],
                            'dn_scale': dn_values[2],
                            # dS (synonymous) parameters
                            'ds_class': ds_values[0],
                            'ds_dist': dist,
                            'ds_shape': ds_values[1],
                            'ds_scale': ds_values[2]
                        }

    return orf_dict

def create_run_dict(orfs_location, kappa, pi, mu, params = None, global_rate = 0.05):
    """
    create dictionary with all keys as required for an HexSE run
    params: for each orf, parameters regarding dn and ds class, shape, dist and scale 

    """
    orfs_info = {}
    for location in orfs_location:
        string_parts = []
        for part in location:
            string_parts.append(','.join([str(value) for value in part]))
            # Include all fragments on orf coordinates separated by colon
            loc_string = ";".join(string_parts)

            # If parameter values for dN and dS are not specified
            if not params:
                dn_values = generate_param_values()
                ds_values = generate_param_values()

            orf_dict = create_orf_dict(loc_string, dn_values, ds_values)
            orfs_info[loc_string] = orf_dict

    info_for_yaml = {
                'global_rate': global_rate,
                'kappa': kappa,
                'pi': pi,
                'mu': mu,
                'orfs': orfs_info
            }

    return info_for_yaml

def main():
    parser = argparse.ArgumentParser(
        description='Create a YAML configuration file to Hexse from GenBank file'
    )

    args = get_args(parser)
    gb_file = args.gb_file
    out = args.yaml_file
    orfs_location = get_locations(gb_file)
    info_for_yaml = create_run_dict(orfs_location, args.k, args.p, args.m)

    with open(out, 'w+') as file:
        yaml.dump(info_for_yaml, file)


if __name__ == '__main__':
    main()
