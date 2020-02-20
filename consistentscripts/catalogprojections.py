#this is for the snapnumber3x2

import yt
import glob
import numpy as np
import random
import os, argparse
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from satellite_analysis.catalogreaders import consistentcatalogreader as consistent 
from satellite_analysis.catalogreaders import tomercatalogreader as tomer

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
                                     description='This is a script for plotting the consistent trees halo catalog onto the 3 2-D coordinate projections of the VELA simulations. This script uses the consistent trees catalog reader found in consistent scripts 
    )
    parser.add_argument('consistent_dir', help='The directory where the consistent trees catalog you want to project is located. The form wanted is like: /nobackupp2/sflarkin/GEN3/rockstar001analysis/VELA07/hlists if you want to run for VELA07')
    parser.add_argument('VELA_dir', help='The directory where the VELA simulation files are located (where the 10MpcBox..., PMcrda0, PMcrs0a0, and stars_a0 files are. The form wanted is like: /nobackupp2/sflarkin/GEN3/VELA07')
    parser.add_argument('VELA_number', help='The VELA number of the catalog you want projected. For VELA 07, you would want to put: 07')
    parser.add_argument('out_dir', help='The output directory where you want the images put. If not current directory, the directory will be created.')
    
    parser.add_argument('--remove_subhalos', nargs='?', default='True', help='Tells the catalog reader whether to remove subhalos or not, default of True')
    parser.add_argument('--halo_mass', nargs='?', default=1e+08, help='Tells the catalog reader the minimum mass of halos in MSun you are interested in. Default of 10**8 MSun.')
    
    args = vars(parser.parse_args())
    return args

args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
VELA_number = args['VELA_number']
out_dir = args['out_dir']
remove_subhalos = args['remove_subhalos']
halo_mass = args['halo_mass']

#check to see if the out_dir is a valid location. if it is not, create it.
if not os.path.exists(out_dir):
    print('Creating output directory')
    os.makedirs(out_dir)
    
#load in the rockstar catatalog halo data
consistent.consistent_catalog_reader(consistent_dir, remove_subhalos=remove_subhalos, halo_mass=halo_mass)

#find the VELA snapshots
VELA_snaps = glob.glob(VELA_dir + '/10MpcBox_csf512_a0*')
VELA_index = []
for snap in VELA_snaps:
    period = [pos for pos, char in enumerate(snap) if char == '.']
    number = snap[period[-2]+1:period[-1]]
    VELA_index.append(number)
VELA_snaps.sort()
VELA_index.sort()

for index in consistent.snapshot_index:
    VELA_a = consistent.consistent_file_index[index]
    print('Generating Graph Grid for snapshot:', VELA_a)
    position = [pos for pos, loc in enumerate(VELA_index) if loc == VELA_a]
    if position == [] or len(position) > 1:
        print('Could not find corresponding VELA 10Mpc File for snapshot:', VELA_a)
    else:
        ds = yt.load(VELA_snaps[position[0]])
        domain_width = float(ds.domain_width.in_units('Mpc/h')[0])
        ad = ds.all_data()
        scale = ds.scale_factor
        masses = yt.np.unique(ad[('darkmatter', 'particle_mass')])
        #filter out the darkmatter0 particles
        def mass_filter(pfilter, data):
            filter = data[(pfilter.filtered_type, 'particle_mass')] == masses[0]
            return filter
        yt.add_particle_filter('darkmatter0', function=mass_filter, filtered_type='darkmatter', requires=['particle_mass'])
        ds.add_particle_filter('darkmatter0')

        x = ad['darkmatter0', "particle_position_x"]
        y = ad['darkmatter0', "particle_position_y"]
        z = ad['darkmatter0', "particle_position_z"]
        
        num_p = int((len(x)/100))
            
        particles_to_plot = random.sample(range(len(x)),int((num_p)))
        
        x0 = [x[index] for index in particles_to_plot]
        y0 = [y[index] for index in particles_to_plot]
        z0 = [z[index] for index in particles_to_plot]

        x1, y1, z1, r1 = [], [], [], [] 
        

        for halos in consistent.halo_data_sorted[index]:
            x1.append(float(halos[17]) * scale / domain_width)
            y1.append(float(halos[18]) * scale / domain_width)
            z1.append(float(halos[19]) * scale / domain_width)
            r1.append(float(halos[11]))
    
    
        print("Successfully generated particle and halo data for snapshot", VELA_a)
        title = '{}/catalogprojection_VELA{}_scale{}.png'.format(out_dir, VELA_number, VELA_a)
    
        plt.figure(figsize=(30,10))
        
        plt.subplot('131')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.scatter(x0, y0, marker='.', s = .1, facecolor='b', edgecolor='b')    
        plt.scatter(x1, y1, s=r1, facecolor='none',edgecolor='g')
        plt.title('Darkmatter + Halos XY Projection')
        plt.legend(('Smallest Darkmatter Particles', 'Consistent Trees Catalog'), markerscale=1)

        #YZ
        plt.subplot('132')
        plt.xlabel('Y')
        plt.ylabel('Z')
        plt.scatter(y0, z0, marker='.', s = .1, facecolor='b', edgecolor='b')    
        plt.scatter(y1, z1, s=r1, facecolor='none',edgecolor='g')
        plt.title('Darkmatter + Halos YZ Projection')

        
        #ZX
        plt.subplot('133')
        plt.xlabel('Z')
        plt.ylabel('X')
        plt.scatter(z0, x0, marker='.', s = .1, facecolor='b', edgecolor='b')    
        plt.scatter(z1, x1, s=r1, facecolor='none',edgecolor='g')
        plt.title('Darkmatter + Halos ZX Projection')

        
        plt.savefig(title)
        print('Saving plot of snap', VELA_a, 'to', title)
        plt.close()
