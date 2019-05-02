#this is for the snapnumber3x2

import yt
import glob
import numpy as np
import random
import argparse
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from satellite_analysis import rockstarcatalogreader as reader 
from satellite_analysis import readtomer as tomer

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('VELA_number')
    args = vars(parser.parse_args())
    return args

args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
VELA_number = args['VELA_number']

#load in the rockstar catatalog halo data
reader.rockstar_catalog_reader(input_dir)

#find the VELA snapshots
VELA_snaps = glob.glob(VELA_dir + '10MpcBox*')
VELA_index = []
for snap in VELA_snaps:
    period = [pos for pos, char in enumerate(snap) if char == '.']
    number = snap[period[-2]+1:period[-1]]
    VELA_index.append(number)
VELA_snaps.sort()
VELA_index.sort()

tomer.read_tomer(VELA_number)
tomer_scales = [str(scale)[2:5] for scale in tomer.tomer_list['aexpn'].tolist()]
r_vir_tomer = tomer.tomer_list['r_vir[kpc]'].tolist()
x_tomer = tomer.tomer_list['center[0](code)'].tolist()
y_tomer = tomer.tomer_list['center[1](code)'].tolist()
z_tomer = tomer.tomer_list['center[2](code)'].tolist()

for index in reader.snapshot_index:
    VELA_a = reader.rockstar_file_index[index]
    print('Generating Graph Grid for snapshot:', VELA_a)
    position = [pos for pos, loc in enumerate(VELA_index) if loc == VELA_a]
    if position == [] or len(position) > 1:
        print('Could not find corresponding VELA 10Mpc File for snapshot:', VELA_a)
    else:
        ds = yt.load(VELA_snaps[position[0]])
        domain_width = float(ds.domain_width.in_units('Mpc/h')[0])
        ad = ds.all_data()
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
            
        particles_to_plot = random.sample(range(len(x)),int((num_p/10)))
        
        x0 = [x[index] for index in particles_to_plot]
        y0 = [y[index] for index in particles_to_plot]
        z0 = [z[index] for index in particles_to_plot]

        #find the corresponding snapshot in tomer's data if it exists
        tomer_number = [pos for pos, number in enumerate(tomer_scales) if number == VELA_a]
        tomer_check = False
        if tomer_number == [] or len(tomer_number) > 1:
            print('Could not find corresponding Tomer Snapshot Data for snap', VELA_a)
            tomer_check = False
        else:
            tomer_check = True
            xt, yto, zt = x_tomer[tomer_number[0]], y_tomer[tomer_number[0]], z_tomer[tomer_number[0]]
            rt = r_vir_tomer[tomer_number[0]]

        #now add the circles for the largest rockstar halos
        #this is the data of the 20 largest halos
        mvir_list_sorted = reader.halo_data_largest[index][1]
        #this is the data for all of the halos
        mvir_list_all = reader.halo_data_all[index]
        x1, x2 = [], []
        y1, y2 = [], []
        z1, z2 = [], []
        r1, r2 = [], []

        for halo0 in mvir_list_sorted:
            x1.append(float(halo0[8])/domain_width)
            y1.append(float(halo0[9])/domain_width)
            z1.append(float(halo0[10])/domain_width)
            r1.append(float(halo0[4]))
    
        for halo1 in mvir_list_all:
            x2.append(float(halo1[8])/domain_width)
            y2.append(float(halo1[9])/domain_width)
            z2.append(float(halo1[10])/domain_width)
            r2.append(float(halo1[4]))
    
        print("Successfully generated particle and halo data for snapshot", VELA_a)
        title = str(VELA_a) +'2x3.png'
    
        plt.figure(figsize=(30,20))
        
        plt.subplot('231')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.scatter(x0, y0, marker='.', s = .1)    
        plt.scatter(x1, y1, s=r1, facecolor='none',edgecolor='tab:orange')
        if tomer_check == True:
            plt.scatter(xt, yto, s=rt, facecolor='none', edgecolor='r')
        plt.title('20 Largest Halos XY Projection')
        plt.legend(('dm0', 'Rockstar Catalog', 'Tomer Center'), markerscale=1)
        
        plt.subplot('234')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.scatter(x0, y0, marker='.', s = .1)    
        plt.scatter(x2, y2, s=r2, facecolor='none',edgecolor='tab:orange')
        plt.title('All Halos XY Projection')
        if tomer_check == True:
            plt.scatter(xt, yto, s=rt, facecolor='none', edgecolor='r')

        #YZ
        plt.subplot('232')
        plt.xlabel('Y')
        plt.ylabel('Z')
        plt.scatter(y0, z0, marker='.', s = .1)    
        plt.scatter(y1, z1, s=r1, facecolor='none',edgecolor='tab:orange')
        plt.title('20 Largest Halos YZ Projection')
        if tomer_check == True:
            plt.scatter(yto, zt, s=rt, facecolor='none', edgecolor='r')

        plt.subplot('235')
        plt.xlabel('Y')
        plt.ylabel('Z')
        plt.scatter(y0, z0, marker='.', s = .1)    
        plt.scatter(y2, z2, s=r2, facecolor='none',edgecolor='tab:orange')
        plt.title('All Halos YZ Projection')
        if tomer_check == True:
            plt.scatter(yto, zt, s=rt, facecolor='none', edgecolor='r')
        
        #ZX
        plt.subplot('233')
        plt.xlabel('Z')
        plt.ylabel('X')
        plt.scatter(z0, x0, marker='.', s = .1)    
        plt.scatter(z1, x1, s=r1, facecolor='none',edgecolor='tab:orange')
        plt.title('20 Largest Halos ZX Projection')
        if tomer_check == True:
            plt.scatter(zt, xt, s=rt, facecolor='none', edgecolor='r')

        plt.subplot('236')
        plt.xlabel('z')
        plt.ylabel('x')
        plt.scatter(z0, x0, marker='.', s = .1)    
        plt.scatter(z2, x2, s=r2, facecolor='none',edgecolor='tab:orange')
        plt.title('All Halos ZX Projection')
        if tomer_check == True:
            plt.scatter(zt, xt, s=rt, facecolor='none', edgecolor='r')
        
        plt.savefig(title)
        print('Saving plot of snap', VELA_a, 'to', title)
        plt.close()
