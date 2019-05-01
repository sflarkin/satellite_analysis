import yt
import glob
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
from mpl_toolkits.axes_grid1 import AxesGrid
from satellite_analysis import rockstarcatalogreader
from satellite_analysis.graphs import halovirplots
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

#need to add arguments for VELA_dir and input_dir for the rockstarcatalogreader
#also add one for making plots of the centers of each halo or for just making a dataset of
#the stellar mass vs darkmatter mass

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--gen_xyz', nargs='?', default=False)
    args = vars(parser.parse_args())
    return args

args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
out_dir = args['out_dir']
gen_xyz = args['--gen_xyz']

#Check if the output directory exists and if it does not, create it.
if not os.path.exists(out_dir):
    print('Creating output directory')
    os.makedirs(out_dir)
if gen_xyz == True:
    xyz_dir = '%s/centergraphs' % outdir
    if not os.path.exists(xyz_dir):
        print('Creating output directory for XYZ graphs')
        os.makedirs(xyz_dir)

#Now run the rockstarcatalogreader to get the halo location data.

#input_dir = '/Users/user1/Documents/VELA07/baserockstar_ascii/subhalos/'
rockstarcatalogreader.rockstar_catalog_reader(input_dir, subhalos=True)

#Collect the names of the VELA simulations in the VELA_dir so that star data can be extracted.
VELA_snaps = glob.glob(VELA_dir + '10MpcBox*')
VELA_snaps.sort()

#Extract the scale factors from the VELA snapshot files to match to the rockstar catalogs.
VELA_index = []
for snap in VELA_snaps:
    period = [pos for pos, char in enumerate(snap) if char == '.']
    number = snap[period[-2]+1:period[-1]]
    VELA_index.append(number)
VELA_index.sort()

#Loop over the rockstar catalogs for each snapshot.
for index in rockstarcatalogreader.snapshot_index:
    VELA_a = rockstarcatalogreader.rockstar_file_index[index]
    position = [pos for pos, loc in enumerate(VELA_index) if loc == VELA_a]
    if position == [] or len(position) > 1:
        print('Could not find corresponding VELA 10Mpc File for snapshot:', VELA_a)
    else:
        print('Finding MVir Masses for snap:', VELA_a)
        ds = yt.load(VELA_snaps[position[0]])
        halo_data = rockstarcatalogreader.halo_data_largest[index][1]
        
        #Create the lists to hold the halo data untill they are written to an ascii file.
        Id_list = []
        stellar_mass_list = []
        darkmatter_mass_list = []
        mass_ratio_list = []
        
        #Loop over the halos to get the center coordinates. Then cut out the spherical region of the simulation
        #with radius of the halos virial radius at the center of the halo, and find the mass of darkmatter
        #and star particles in that region.
        for halo in halo_data:
            domain_width = float(ds.domain_width.in_units('Mpc/h')[0])
            print(domain_width)
            x = float(halo[8])/domain_width
            y = float(halo[9])/domain_width
            z = float(halo[10])/domain_width
            rvir = float(halo[4])
            Id = halo[0]
            center = [x, y, z]
            sp = ds.sphere(center, (rvir, 'kpc/h'))
            stellar_mass, darkmatter_mass = sp.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                          ('darkmatter', 'particle_mass')])
            
            #Add the halo values to the lists to be written to an ascii file.
            Id_list.appned(Id)
            stellar_mass_list.append(stallar_mass)
            darkmatter_mass_list.append(darkmatter_mass)
            mass_ratio_list.append(stellar_mass/darkmatter_mass)
            
            #Now make the xyz graphs if wanted.
            if gen_xyz == True:
                #check if the current Id is in the top X halos we want graphs for
                halo_vir_plots from graphs folder
            
        #Now write the halo mass information to an ascii file
        file_name = 
        data = Table([Id_list, stellar_mass_list, darkmatter_mass_list, mass_ratio_list],\
                     names=['Id', 'stellar_mass', 'darkmatter_mass', 'mass_ratio'])
        ascii.write(data, output=file_name)
            
