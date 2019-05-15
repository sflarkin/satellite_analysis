import yt
import numpy as np
import argparse
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import glob

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('VELA_dir')
    args = vars(parser.parse_args())
    return args

args = parse()

VELA_dir = args['VELA_dir']

all_files = glob.glob('%s/10Mpc*' % VELA_dir)
all_files.sort()
domain_width, domain_dimensions = [], []
smallest_dx, max_level, current_redshift = [], [], []

for file in all_files:
    ds = yt.load(file)
    domain_width.append(ds.domain_width.in_units('Mpc/h')[0])
    domain_dimensions.append(ds.domain_dimensions[0])
    max_level.append(ds.max_level)
    smallest_dx.append(ds.index.get_smallest_dx().in_units('Mpc/h'))
    current_redshift.append(ds.current_redshift)

file_name = '%s/smallestgrid.ascii' % VELA_dir
data = Table([current_redshift, smallest_dx, domain_width, domain_dimensions, max_level],\
              names=['current_redshift', 'smallest_dx(pc/h)', 'domain_width(Mpc/h)', 'domain_dimension', 'max_level'])
ascii.write(data, output=file_name)