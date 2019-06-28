import yt
import glob
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
from mpl_toolkits.axes_grid1 import AxesGrid
from satellite_analysis.catalogreaders import consistentcatalogreader as consistent
from satellite_analysis.graphs import stellarmassrelationfindingrvirprojections as rvirprj
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
    parser.add_argument('--gen_xyz', nargs='?', default='False')
    parser.add_argument('--subhalos', nargs='?', default='False')
    args = vars(parser.parse_args())
    return args

args = parse()
print(args)

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
out_dir = args['out_dir']
gen_xyz = args['gen_xyz']
subhalos = args['subhalos']

if gen_xyz == 'True':
    print('Generating Images')
#Check if the output directory exists and if it does not, create it.
if not os.path.exists(out_dir):
    print('Creating output directory')
    os.makedirs(out_dir)
#if gen_xyz == True:
#    xyz_dir = '%s/centergraphs' % outdir
#    if not os.path.exists(xyz_dir):
#        print('Creating output directory for XYZ graphs')
#        os.makedirs(xyz_dir)

#Now run the rockstarcatalogreader to get the halo location data.

#input_dir = '/Users/user1/Documents/VELA07/baserockstar_ascii/subhalos/'
consistent.consistent_catalog_reader(input_dir, subhalos=subhalos, halo_number=1)

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

#Loop over the rockstar catalogs for each snapshot.
for index in consistent.snapshot_index:
    VELA_a = consistent.consistent_file_index[index]
    position = [pos for pos, loc in enumerate(VELA_index) if loc == VELA_a]
    if position == [] or len(position) > 1:
        print('Could not find corresponding VELA 10Mpc File for snapshot:', VELA_a)
    else:
        print('Finding MVir Masses for snap:', VELA_a)
        ds = yt.load(VELA_snaps[position[0]])
        halo_data = consistent.halo_data_all[index]
        
        #Create the lists to hold the halo data untill they are written to an ascii file.
        Id_list = []
        pid_list = []
        mvir_list = []
        mpeak_list = []
        stellar_mass_10rvir_list = []
        darkmatter_mass_10vir_list = []
        stellar_mass_15rvir_list = []
        darkmatter_mass_15rvir_list = []
        stellar_mass_20vir_list = []
        darkmatter_mass_20rvir_list = []
        stellar_mass_rvir_list = []
        darkmatter_mass_rvir_list = []
        
        #create the list of halo_ids of the largest halos to plot if desired
        if gen_xyz == 'True':
            largest_halo_ids = [consistent.halo_data_largest[index][0][q][1] for q in range(len(consistent.halo_data_largest[index][0]))]
            print(largest_halo_ids)
        #Loop over the halos to get the center coordinates. Then cut out the spherical region of the simulation
        #with radius of the halos virial radius at the center of the halo, and find the mass of darkmatter
        #and star particles in that region.
        for halo in halo_data:
            domain_width = float(ds.domain_width.in_units('Mpc/h')[0])
            x = float(halo[17])/domain_width
            y = float(halo[18])/domain_width
            z = float(halo[19])/domain_width
            rvir = float(halo[11])
            Id = halo[1]
            center = [x, y, z]
            sp = ds.sphere(center, (rvir, 'kpc/h'))
            stellar_mass_rvir, darkmatter_mass_rvir = sp.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                          ('darkmatter', 'particle_mass')])
            sp1 = ds.sphere(center, (rvir*.1, 'kpc/h'))
            stellar_mass_10rvir, darkmatter_mass_10rvir = sp1.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                          ('darkmatter', 'particle_mass')])
            
            sp2 = ds.sphere(center, (rvir*.15, 'kpc/h'))
            stellar_mass_15rvir, darkmatter_mass_15rvir = sp2.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                          ('darkmatter', 'particle_mass')])
            
            sp3 = ds.sphere(center, (rvir*.2, 'kpc/h'))
            stellar_mass_20rvir, darkmatter_mass_20rvir = sp3.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                          ('darkmatter', 'particle_mass')])
            #get the Mpeak, Mvir, and pid if its a satellite or not
            pid = halo[5]
            mvir = halo[10]
            mpeak = halo[61]
            #Add the halo values to the lists to be written to an ascii file.
            Id_list.append(Id)
            pid_list.append(pid)
            mvir_list.append(mvir)
            mpeak_list.append(mpeak)
            stellar_mass_10rvir_list.append(stellar_mass_10rvir.in_units('Msun/h'))
            darkmatter_mass_10vir_list.append(darkmatter_mass_10rvir.in_units('Msun/h'))
            stellar_mass_15rvir_list.append(stellar_mass_15rvir.in_units('Msun/h'))
            darkmatter_mass_15rvir_list.append(darkmatter_mass_15rvir.in_units('Msun/h'))
            stellar_mass_20vir_list.append(stellar_mass_20rvir.in_units('Msun/h'))
            darkmatter_mass_20rvir_list.append(darkmatter_mass_20rvir.in_units('Msun/h'))
            stellar_mass_rvir_list.append(stellar_mass_rvir.in_units('Msun/h'))
            darkmatter_mass_rvir_list.append(darkmatter_mass_rvir.in_units('Msun/h'))
            
            #Now make the xyz graphs if wanted.
            if gen_xyz == 'True':
                #check if the current Id is in the top X halos we want graphs for
                if Id in largest_halo_ids:
                    rvirprj.rvirprojections(ds, center, rvir, sp, out_dir, VELA_a, Id)
    
        #Now write the halo mass information to an ascii file
        file_name = '%s/halomass%s.ascii' % (out_dir, VELA_a)
        data = Table([Id_list, pid_list, mvir_list, mpeak_list, stellar_mass_10rvir_list, darkmatter_mass_10rvir_list, \
                      stellar_mass_15rvir_list, darkmatter_mass_15rvir_list, stellar_mass_20vir_list, \
                      darkmatter_mass_20rvir_list, stellar_mass_rvir_list, darkmatter_mass_rvir_list],\
                     names=['Id[1]', 'Pid[5]', 'Mvir[11]', 'Mpeak[61]', 'stellarmass_.1rvir', 'darkmatter_mass_.1rivr',\
                            'stellarmass_.15rvir', 'darkmattermass_.15rvir', 'stellarmass_.2rvir', 'darkmattermass_.2rvir',\
                           'stellarmass_rvir', 'darkmattermass_rvir'])
        ascii.write(data, output=file_name, overwrite=True)
            
