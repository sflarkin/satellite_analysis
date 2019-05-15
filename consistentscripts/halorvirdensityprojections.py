import yt
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
from mpl_toolkits.axes_grid1 import AxesGrid
from satellite_analysis.catalogreaders import consistentcatalogreader as consistent
from satellite_analysis.catalogreaders import tomercatalogreader as tomer
from satellite_analysis.graphs import stellarmassrelationfindingrvirprojections as rvirprj


#need to add arguments for VELA_dir and input_dir for the rockstarcatalogreader

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('VELA_number')
    parser.add_argument('--subhalos', nargs='?', default=True)
    args = vars(parser.parse_args())
    return args

args = parse()
#print(args)

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
VELA_number = args['VELA_number']
subhalos = args['subhalos']

#input_dir = '/Users/user1/documents/VELA07/hlists'
#VELA_dir = '/Users/user1/documents/VELA07'
#subhalos = False

#Now run the rockstarcatalogreader to get the halo location data.

#input_dir = '/Users/user1/Documents/VELA07/baserockstar_ascii/subhalos/'
consistent.consistent_catalog_reader(input_dir, subhalos=subhalos)

#Collect the names of the VELA simulations in the VELA_dir so that star data can be extracted.
VELA_snaps = glob.glob(VELA_dir + '/10MpcBox*')
VELA_snaps.sort()

#Extract the scale factors from the VELA snapshot files to match to the rockstar catalogs.
VELA_index = []
for snap in VELA_snaps:
    period = [pos for pos, char in enumerate(snap) if char == '.']
    number = snap[period[-2]+1:period[-1]]
    VELA_index.append(number)
VELA_index.sort()

#load the tomer catalogs
tomer.read_tomer(VELA_number)
tomer_scales = [str(scale)[2:5] for scale in tomer.tomer_list['aexpn'].tolist()]
r_vir_tomer = tomer.tomer_list['r_vir[kpc]'].tolist()
x_tomer = tomer.tomer_list['center[0](code)'].tolist()
y_tomer = tomer.tomer_list['center[1](code)'].tolist()
z_tomer = tomer.tomer_list['center[2](code)'].tolist()

#Loop over the rockstar catalogs for each snapshot.
for index in consistent.snapshot_index:
    VELA_a = consistent.consistent_file_index[index]
    position = [pos for pos, loc in enumerate(VELA_index) if loc == VELA_a]
    if position == [] or len(position) > 1:
        print('Could not find corresponding VELA 10Mpc File for snapshot:', VELA_a)
    else:
        print('Finding MVir Masses for snap:', VELA_a)
        #load the yt snap and the halo_data from the catalog reader
        ds = yt.load(VELA_snaps[position[0]])
        halo_data = consistent.halo_data_largest[index][0]
        
        #define the parameters for the largest halo for plotting
        domain_width = float(ds.domain_width.in_units('Mpc/h')[0])
        largest_halo = halo_data[0]
        x = float(largest_halo[17])/domain_width
        y = float(largest_halo[18])/domain_width
        z = float(largest_halo[19])/domain_width
        rvir = float(largest_halo[11])
        center = [x, y, z]
        print(rvir, center)
        #extract the rvir sphere around the largest halo 
        sp = ds.sphere(center, (rvir, 'kpc/h'))
        
        tomer_number = [pos for pos, number in enumerate(tomer_scales) if number == VELA_a]
        if tomer_number == [] or len(tomer_number) > 1:
            print('Could not find corresponding Tomer Snapshot Data for snap', VELA_a)
            tomer_center = [0,0,0]
            tomer_rvir = [0]
        else:
            tomer_center = [x_tomer[tomer_number[0]], y_tomer[tomer_number[0]], z_tomer[tomer_number[0]]]
            tomer_rvir = r_vir_tomer[tomer_number[0]]
        
        #call the plotting function
        rvirprj.rvirprojectionsallhalos(ds, center, rvir, sp, tomer_center, tomer_rvir, halo_data, domain_width, input_dir, VELA_a)
    
            
