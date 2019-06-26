import yt
import numpy as np
import argparse
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import glob

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('VELA_dir')
    parser.add_argument('out_dir')
    args = vars(parser.parse_args())
    return args

args = parse()

VELA_dir = args['VELA_dir']
out_dir = args['out_dir']
all_files = glob.glob(VELA_dir + '/10Mpc*')
completed = glob.glob('*.ascii')
print('Files already converted:', completed)
digits_list = []
for file in all_files:

    #extract the snapshot number for use
    digits = file[-5:-2]
    digits_list.append(str(digits))
    file_name = VELA_dir + '/' + str(digits) + '.ascii'
    print(file_name)
    if file_name not in completed:

        #load the particles for writing

        ds = yt.load(file)

        ad = ds.all_data()
        masses = yt.np.unique(ad[('darkmatter', 'particle_mass')])
        #filter out the darkmatter0 particles
        def mass_filter(pfilter, data):
            filter = data[(pfilter.filtered_type, 'particle_mass')] == masses[0]
            return filter
        yt.add_particle_filter('darkmatter0', function=mass_filter, filtered_type='darkmatter', requires=['particle_mass'])
        ds.add_particle_filter('darkmatter0')

        x = ad['darkmatter0', "particle_position_x"].in_units('Mpc/h')
        y = ad['darkmatter0', "particle_position_y"].in_units('Mpc/h')
        z = ad['darkmatter0', "particle_position_z"].in_units('Mpc/h')

        vx = ad['darkmatter0', "particle_velocity_x"].in_units('km/s')
        vy = ad['darkmatter0', "particle_velocity_y"].in_units('km/s')
        vz = ad['darkmatter0', "particle_velocity_z"].in_units('km/s')

        ids = [q for q in range(len(x))]

        #create the scale factor line for the ascii comment section
        scale_factor = str(ds.scale_factor)
        a_line = 'a = ' + scale_factor
        comment = {'comments':[a_line]}

        #create the ascii table
        data = Table([x, y, z, vx, vy, vz, ids], names=['x', 'y', 'z', 'vx', 'vy', 'vz', 'id'], meta=comment)

        ascii.write(data, output=file_name, comment='#')

#now write the snapshot_names.txt file

digits_list.sort()

f = open('%s/snapshot_names.txt' % out_dir, 'w')
for item in digits_list:
    f.write("%s\n" % item)
f.close()
