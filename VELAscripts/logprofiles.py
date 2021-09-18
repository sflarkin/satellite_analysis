import yt
import numpy as np
from yt import YTQuantity

import glob
import os, argparse

from satellite_analysis.catalogreaders import consistentcatalogreader as consistent

from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('proper_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--remove_subhalos', nargs='?', default='True')
    parser.add_argument('--halo_mass', nargs='?', default=10e+09)
    args = vars(parser.parse_args())
    return args


args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
proper_dir = args['proper_dir']
out_dir = args['out_dir']
VELA_snaps = glob.glob(VELA_dir + '/10MpcBox*')
VELA_snaps.sort()
print(VELA_snaps)

consistent.consistent_catalog_reader(input_dir, remove_subhalos='True', halo_mass=10e9)

completed = glog.glob('{}/logspacedprofiles_scale*_id*.ascii'.format(out_dir))

for snap in reversed(VELA_snaps):
    ds = yt.load(snap)
    current_scale = str(ds.scale_factor)[2:5]
    
    profile_file_name ='{}/logspacedprofiles_scale{}_id{}.ascii'.format(out_dir, current_scale, rockstar_id)
    
    if profile_file_name in completed:
        continue
    print('Trying to find matching Consistent File for scale: {}'.format(current_scale))
    
    position = [pos for pos, loc in enumerate(consistent.consistent_file_index) if loc == current_scale]

    if position == [] or len(position) > 1:
        print('Could not find corresponding consistent trees data for file {}'.format(snap))
    else:
        
        print('Finding MVir Masses for snap:{}'.format(snap))
        print('Matching Scales Test:{} {}'.format(current_scale, consistent.consistent_file_index[position[0]]))
        
        #make the lists to hold the data to be written to files later
        count = 0
        
        
        halo_data = consistent.halo_data_sorted[position[0]]

        for halo in halo_data:
            
            #set the values to be calculated to None, so that if any values arent found we will know
            
            properfile = '{}/improvedrvir{}.ascii'.format(proper_dir, current_scale)
            
            
            rockstar_id = float(halo[1])
            print(rockstar_id)
            x = ds.quan(float(halo[17]), 'Mpccm/h')
            y = ds.quan(float(halo[18]), 'Mpccm/h')
            z = ds.quan(float(halo[19]), 'Mpccm/h')
            center = [x,y,z]
            
            #now we need to read the rvir data 
            
            readfile = open(properfile)
            lines = readfile.readlines()
            
            improved_rvir = 0
            
            for line in lines:
                #ignore the first line, as it is only text explaing what each column is
                if line == lines[0]:
                    continue
                    
                current_line = line.split()
                #print(current_line[0])
                if float(current_line[0]) == rockstar_id:
                
                    improved_rvir = current_line[3]
                
            rvir = ds.quan(float(improved_rvir), 'kpc')
            print(rvir)

            
            sp = ds.sphere(center, rvir)
            space = np.append(np.array([0]), np.logspace(-2, 0, 30))

            radius_bins = ds.arr(space*float(rvir), 'kpc') 
            #radius_bins = ds.arr(np.logspace(-2, 0, 30)*float(rvir), 'kpc')
            #print(radius_bins)
            
            #add set varialbles to none so if they do not get set, we see that when the get written
                     

            rp_gas = yt.create_profile(sp, 'radius', [('gas', 'cell_mass')], accumulation = True, 
                                        units = {'radius': 'kpc', 'cell_mass': 'Msun'}, weight_field=None,
                                        override_bins={'radius': radius_bins})

            rp_stars = yt.create_profile(sp, ('stars', 'particle_radius'), [('stars', 'particle_mass')], accumulation = True,
                                         units = {('stars', 'particle_radius'): 'kpc', ('stars', 'particle_mass'): 'Msun'},
                                         weight_field=None, override_bins={('stars', 'particle_radius'): radius_bins})
            
            rp_darkmatter = yt.create_profile(sp, ('darkmatter', 'particle_radius'), [('darkmatter', 'particle_mass')], accumulation = True,
                                              units = {('darkmatter', 'particle_radius'): 'kpc', ('darkmatter', 'particle_mass'): 'Msun'},
                                              weight_field=None, override_bins={('darkmatter', 'particle_radius'): radius_bins})
            
            
            gas_profile_list = []
            star_profile_list = []
            darkmatter_profile_list = []
            profile_radii_list = []
            
            for bins in range(len(rp_gas[('gas', 'cell_mass')])):
                gas_radius = rp_gas.x_bins[bins+1]
                gas_mass = rp_gas[('gas', 'cell_mass')][bins]
    
                star_radius = rp_stars.x_bins[bins+1]
                star_mass = rp_stars[('stars', 'particle_mass')][bins]
    
                darkmatter_radius = rp_darkmatter.x_bins[bins+1]
                darkmatter_mass = rp_darkmatter[('darkmatter', 'particle_mass')][bins]
                
                #write the profile to lists for writing to a file
                gas_profile_list.append(gas_mass.in_units('Msun'))
                star_profile_list.append(star_mass.in_units('Msun'))
                darkmatter_profile_list.append(darkmatter_mass.in_units('Msun'))
                profile_radii_list.append(gas_radius.in_units('kpc'))
                

                
            #write the lists to a file
            profile_file_name ='{}/logspacedprofiles_scale{}_id{}.ascii'.format(out_dir, current_scale, rockstar_id)
            profile_data = [profile_radii_list, gas_profile_list, star_profile_list, darkmatter_profile_list]
            profile_names = ['Radius_kpc', 'Gas_Mass_Msun', 'Stellar_Mass_Msun', 'Darkmatter_Mass_Msun']
            profile_data = Table(profile_data, names=profile_names)
            ascii.write(profile_data, output=profile_file_name, overwrite=True)
