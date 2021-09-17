
import to_yaml

def get_settings(args):
    """
    """
    yaml_settings = to_yaml([args.seq, csv, args.config])
    settings = {
       **vars(args)
       yaml_settings
    }
    print (settings)
    return settings