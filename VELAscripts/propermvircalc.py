import yt
import numpy as np
from yt import YTQuantity

import glob
import os, argparse
from yt.utilities.cosmology import Cosmology

from satellite_analysis.catalogreaders import consistentcatalogreader as consistent
from satellite_analysis.graphs import stellarmassrelationfindingrvirprojections as rvirprj

from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--remove_subhalos', nargs='?', default='True')
    parser.add_argument('--halo_mass', nargs='?', default=10e+09)
    args = vars(parser.parse_args())
    return args

def delta_vir_calculator(z, omega_matter, omega_lambda):
    OmegaA= omega_matter / ( omega_matter + omega_lambda* (1/(1+z))**3)
    deltavir = (18*np.pi**2 + 82*(OmegaA -1) - 39*(OmegaA-1)**2 )/OmegaA
    return deltavir

args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
out_dir = args['out_dir']
remove_subhalos = args['remove_subhalos']
halo_mass = args['halo_mass']



#Check if the output directory exists and if it does not, create it.
if not os.path.exists(out_dir):
    print('Creating output directory')
    os.makedirs(out_dir)

VELA_snaps = glob.glob(VELA_dir + '/10MpcBox*')
VELA_snaps.sort()
print(VELA_snaps)

consistent.consistent_catalog_reader(input_dir, remove_subhalos='True', halo_mass=10e9)

completed = glob.glob('{}/improvedrvir*.ascii'.format(out_dir))

for snap in reversed(VELA_snaps):
    ds = yt.load(snap)
    current_scale = str(ds.scale_factor)[2:5]
    

    file_name = '{}/improvedrvir{}.ascii'.format(out_dir, current_scale)
    
    if file_name in completed:
        print('Already Analyzed Snap')
        continue
              
        
    print('Trying to find matching Consistent File for scale: {}'.format(current_scale))
    
    position = [pos for pos, loc in enumerate(consistent.consistent_file_index) if loc == current_scale]

    if position == [] or len(position) > 1:
        print('Could not find corresponding consistent trees data for file {}'.format(snap))
    else:
        
        print('Finding MVir Masses for snap:{}'.format(snap))
        print('Matching Scales Test:{} {}'.format(current_scale, consistent.consistent_file_index[position[0]]))
        
        #make the lists to hold the data to be written to files later
        
        rockstar_rvir_list = []
        rockstar_mvir_list = []
        rockstar_id_list = []
        
        new_rvir_list = []
        new_rvir_delta_list = []
        
        stellar_mass_rvir_list, gas_mass_rvir_list, darkmatter_mass_rvir_list = [], [], []
        stellar_mass_010_list, gas_mass_010_list = [], []
        stellar_mass_015_list, gas_mass_015_list = [], []
        stellar_mass_020_list, gas_mass_020_list = [], []
        
        stellar_halfmass_010_list, stellar_halfmass_015_list, stellar_halfmass_020_list = [], [], []
        stellar_halfmass_010_delta_list, stellar_halfmass_015_delta_list, stellar_halfmass_020_delta_list = [], [], []
        
        count = 0
        
        
        halo_data = consistent.halo_data_sorted[position[0]]
        print(halo_data)

        for halo in halo_data:
            
            #set the values to be calculated to None, so that if any values arent found we will know
            rvir, mvir, rockstar_id = None, None, None
            new_rvir, new_rvir_delta = None, None
            stellar_mass_rvir, gas_mass_rvir, darkmatter_mass_rvir = None, None, None
            stellar_mass_010, stellar_mass_015, stellar_mass_020 = None, None, None
            gas_mass_010, gas_mass_015, gas_mass_020 = None, None, None
            stellar_halfmass_010, stellar_halfmass_015, stellar_halfmass_020 = None, None, None
            stellar_halfmass_010_delta, stellar_halfmass_015_delta, stellar_halfmass_020_delta = None, None, None
            
            
            
            x = ds.quan(float(halo[17]), 'Mpccm/h')
            y = ds.quan(float(halo[18]), 'Mpccm/h')
            z = ds.quan(float(halo[19]), 'Mpccm/h')
            mvir = ds.quan(float(halo[10]), 'Msun/h')
            rvir = ds.quan(float(halo[11]), 'kpccm/h')
            rockstar_id = float(halo[1])
            center = [x, y, z]
            print(rvir.in_units('kpc'), center)
            
            sp_proj = ds.sphere(center, (float(rvir.in_units('kpc')), 'kpc'))
            
            rvirprj.rvirprojections(ds, center, rvir, sp_proj, out_dir, current_scale, rockstar_id, count)
            count = count+1
            
            sp = ds.sphere(center, (float(rvir.in_units('kpc'))*2.25, 'kpc'))
            
            radius_bins = ds.arr(np.linspace(0, float(rvir.in_units('kpc'))*2.25, num = 100), 'kpc')
            print(rockstar_id)
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
            
            
        
            #now loop through the bins to find the right delta_vir
            z = ds.current_redshift
            a = 1 / (1+z)
            omega_matter = ds.omega_matter
            omega_lambda = ds.omega_lambda
            h_0 = YTQuantity(0.699999988079071, '100*km/s/Mpc')
            G1 = YTQuantity(6.67430e-11, 'm*m*m/kg/s/s')
            
            delta_vir_actual = delta_vir_calculator(z, omega_matter, omega_lambda)
            print('Wanted detla_vir: {}'.format(delta_vir_actual))
            #now use the omega_matter_current, and the critical density to find the current matter_density
            critical_density = 3*((h_0)**2) / (8 * np.pi * G1) * a**-3

            rho_matter = omega_matter * critical_density
        
            
            #this is where all of the differences of delta vir will be stored so that the clocest can be saved and recalculated
            delta_vir_distance = []
            
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
    
    
                total_mass = gas_mass + star_mass + darkmatter_mass
    
                delta_vir_current = (total_mass / gas_radius**3) * (3 / (4 * np.pi)) * (1/rho_matter.in_units('Msun/kpc**3'))
                print(gas_radius, total_mass, delta_vir_current)
                delta_vir_distance.append(float(delta_vir_current - delta_vir_actual))
                
                #write the profile to lists for writing to a file
                gas_profile_list.append(gas_mass.in_units('Msun'))
                star_profile_list.append(star_mass.in_units('Msun'))
                darkmatter_profile_list.append(darkmatter_mass.in_units('Msun'))
                profile_radii_list.append(gas_radius.in_units('kpc'))
                

                
            #write the lists to a file
            profile_file_name ='{}/improvedrvirprofiles_scale{}_id{}.ascii'.format(out_dir, current_scale, rockstar_id)
            profile_data = [profile_radii_list, gas_profile_list, star_profile_list, darkmatter_profile_list]
            profile_names = ['Radius_kpc', 'Gas_Mass_Msun', 'Stellar_Mass_Msun', 'Darkmatter_Mass_Msun']
            profile_data = Table(profile_data, names=profile_names)
            ascii.write(profile_data, output=profile_file_name, overwrite=True)
            
            
            print(delta_vir_distance)
            didnot = False
            for index in range(len(delta_vir_distance)):
                if index == range(len(delta_vir_distance))[-1]:
        
                    print('Did not find Delta Vir Crossover point')
                    didnot = True
                else:
                    if delta_vir_distance[index] * delta_vir_distance[index+1] < 0:
                        #when this is true, we have the right sphere limits to expand upon
                        gas_to_add = rp_gas[('gas', 'cell_mass')][index]
                        stars_to_add = rp_stars[('stars', 'particle_mass')][index]
                        darkmatter_to_add = rp_darkmatter[('darkmatter', 'particle_mass')][index]
                        radius1 = rp_gas.x_bins[index+1]
                        radius2 = rp_gas.x_bins[index+2]
            
                        radius_bins_refined_delta_vir = ds.arr(np.linspace(radius1, radius2, num = 100), 'kpc')
            
            
                        rp_gas_refined = yt.create_profile(sp, 'radius', [('gas', 'cell_mass')], accumulation = True, 
                                                    units = {'radius': 'kpc', 'cell_mass': 'Msun'}, weight_field=None,
                                                    override_bins={'radius': radius_bins_refined_delta_vir})

                        rp_stars_refined = yt.create_profile(sp, ('stars', 'particle_radius'), [('stars', 'particle_mass')], accumulation = True,
                                                     units = {('stars', 'particle_radius'): 'kpc', ('stars', 'particle_mass'): 'Msun'},
                                                     weight_field=None, override_bins={('stars', 'particle_radius'): radius_bins_refined_delta_vir})
            
                        rp_darkmatter_refined = yt.create_profile(sp, ('darkmatter', 'particle_radius'), [('darkmatter', 'particle_mass')], accumulation = True,
                                                          units = {('darkmatter', 'particle_radius'): 'kpc', ('darkmatter', 'particle_mass'): 'Msun'},
                                                          weight_field=None, override_bins={('darkmatter', 'particle_radius'): radius_bins_refined_delta_vir})
           
                        delta_vir_distance_refined = []
            
                        for refined_bins in range(len(rp_gas_refined[('gas', 'cell_mass')])):
                    
                    
                            if refined_bins == 0:
                                gas_radius_refined = rp_gas_refined.x_bins[refined_bins]
                                total_mass_refined = gas_to_add + stars_to_add + darkmatter_to_add
                                delta_vir_refined = (total_mass_refined / gas_radius_refined**3) * (3 / (4 * np.pi)) * (1/rho_matter.in_units('Msun/kpc**3'))
                                delta_vir_distance_refined.append(float(delta_vir_refined - delta_vir_actual))
                                print(gas_radius_refined, total_mass_refined, delta_vir_refined)
                            
                                
                            gas_radius_refined = rp_gas_refined.x_bins[refined_bins+1]
                            gas_mass_refined = rp_gas_refined[('gas', 'cell_mass')][refined_bins] + gas_to_add
    
                            star_radius_refined = rp_stars_refined.x_bins[refined_bins+1]
                            star_mass_refined = rp_stars_refined[('stars', 'particle_mass')][refined_bins] + stars_to_add
    
                            darkmatter_radius_refined = rp_darkmatter_refined.x_bins[refined_bins+1]
                            darkmatter_mass_refined = rp_darkmatter_refined[('darkmatter', 'particle_mass')][refined_bins] + darkmatter_to_add
    
                            total_mass_refined = gas_mass_refined + star_mass_refined + darkmatter_mass_refined
    
                            delta_vir_refined = (total_mass_refined / gas_radius_refined**3) * (3 / (4 * np.pi)) * (1/rho_matter.in_units('Msun/kpc**3'))
                            delta_vir_distance_refined.append(float(delta_vir_refined - delta_vir_actual))
                            print(gas_radius_refined, total_mass_refined, delta_vir_refined)
                        #print(delta_vir_distance_refined)
            
                        for index_refined in range(len(delta_vir_distance_refined)):
                            if index_refined == range(len(delta_vir_distance_refined))[-1]:
                                
                                print('Did not find Delta Vir Crossover point')
                                break
                            else:
                                if delta_vir_distance_refined[index_refined] * delta_vir_distance_refined[index_refined+1] < 0:
                                    radius1_refined = rp_gas_refined.x_bins[index_refined+0]
                                    radius2_refined = rp_gas_refined.x_bins[index_refined+1]
                                    new_rvir = (radius1_refined + radius2_refined) / 2
                                    new_rvir_delta = new_rvir - radius1_refined
                                    print(radius1_refined, radius2_refined)
                                    break
                        break
            
            if didnot == True:
                continue
            
            #now that we have the new correct Mvir, we collect the values we want to store for it 
        
            #first, find the stars and gas within .1, .15, and .2 Rvir, and the darkmatter withing 1 Rvir
            sp_new_rvir = ds.sphere(center, (float(new_rvir.in_units('kpc')), 'kpc'))
            sp_new_rvir_010 = ds.sphere(center, (float(new_rvir.in_units('kpc'))*.1, 'kpc'))
            sp_new_rvir_015 = ds.sphere(center, (float(new_rvir.in_units('kpc'))*.15, 'kpc'))
            sp_new_rvir_020 = ds.sphere(center, (float(new_rvir.in_units('kpc'))*.2, 'kpc'))
            
            
            
            stellar_mass_rvir, darkmatter_mass_rvir, gas_mass_rvir = sp_new_rvir.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                                               ('darkmatter', 'particle_mass'),\
                                                                                               ('gas', 'cell_mass')])
            stellar_mass_010, gas_mass_010 = sp_new_rvir_010.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                                               ('gas', 'cell_mass')])

            stellar_mass_015, gas_mass_015 = sp_new_rvir_015.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                                               ('gas', 'cell_mass')])
            
            stellar_mass_020, gas_mass_020 = sp_new_rvir_020.quantities.total_quantity([('stars', 'particle_mass'),\
                                                                                               ('gas', 'cell_mass')])
            
            #now find the half_mass radius for the stars and gas of each sub radius
            
            stellar_halfmass_010_difference, stellar_halfmass_015_difference, stellar_halfmass_020_difference = [], [], []
            gas_halfmass_010_difference, gas_halfmass_015_difference, gas_halfmass_020_difference = [], [], []
            
            for index_3 in range(len(rp_stars[('stars', 'particle_mass')])):
                
                if index_3 == 0:
                    stellar_halfmass_010_difference.append(stellar_mass_010/2-0)
                    stellar_halfmass_015_difference.append(stellar_mass_015/2)
                    stellar_halfmass_020_difference.append(stellar_mass_020/2)
                    
                stellar_halfmass_010_difference.append(stellar_mass_010/2 - rp_stars[('stars', 'particle_mass')][index_3])
                stellar_halfmass_015_difference.append(stellar_mass_015/2 - rp_stars[('stars', 'particle_mass')][index_3])
                stellar_halfmass_020_difference.append(stellar_mass_020/2 - rp_stars[('stars', 'particle_mass')][index_3])
                
                gas_halfmass_010_difference.append(gas_mass_010/2 - (rp_gas[('gas', 'cell_mass')][index_3]))
                gas_halfmass_015_difference.append(gas_mass_015/2 - (rp_gas[('gas', 'cell_mass')][index_3]))
                gas_halfmass_020_difference.append(gas_mass_020/2 - (rp_gas[('gas', 'cell_mass')][index_3]))
        
            
            print('Stellar_Halfmass_010_Difference')
            print(stellar_halfmass_010_difference)
            print('Stellar_Halfmass_015_Difference')
            print(stellar_halfmass_015_difference)
            print('Stellar_Halfmass_020_Difference')
            print(stellar_halfmass_020_difference)
            
            print('Wanted Stellar010: {}'.format(stellar_mass_010.in_units('g')/2))
            print('Wanted Stellar015: {}'.format(stellar_mass_015.in_units('g')/2))
            print('Wanted Stellar020: {}'.format(stellar_mass_020.in_units('g')/2))
            
            for index_4 in range(len(stellar_halfmass_010_difference)):
                print(index_4, stellar_halfmass_010_difference[index_4], stellar_halfmass_015_difference[index_4], stellar_halfmass_020_difference[index_4])
                if index_4 == range(len(stellar_halfmass_010_difference))[-1]:
        
                    print('Did not find Halfmass Crossover point')
                    break
                else:
        
                    if stellar_halfmass_010_difference[index_4] * stellar_halfmass_010_difference[index_4+1] < 0:
                        print('Finding 010 halfmass')
                        if index_4 == 0:
                            stars_halfmass_to_add = 0
                            print('0')
                        else:
                            stars_halfmass_to_add = rp_stars[('stars', 'particle_mass')][index_4-1]
                            print(stars_halfmass_to_add.in_units('g'))
                        
                        radius_3 = rp_stars.x_bins[index_4+0]
                        radius_4 = rp_stars.x_bins[index_4+1]
                        
                        radius_bins_refined_halfmass = ds.arr(np.linspace(radius_3, radius_4, num = 100), 'kpc')

                        rp_stars_halfmass = yt.create_profile(sp, ('stars', 'particle_radius'), [('stars', 'particle_mass')], accumulation = True,
                                                                 units = {('stars', 'particle_radius'): 'kpc', ('stars', 'particle_mass'): 'Msun'},
                                                                 weight_field=None, override_bins={('stars', 'particle_radius'): radius_bins_refined_halfmass})

                        stellar_halfmass_010_difference_refined = []

                        for halfmass_refined_bins in range(len(rp_stars_halfmass[('stars', 'particle_mass')])):
                            #star_radius_refined = rp_stars_halfmass.x_bins[halfmass_refined_bins + 1]
                            if halfmass_refined_bins == 0:
                                stellar_halfmass_010_difference_refined.append(stellar_mass_010/2 - stars_halfmass_to_add)
                                
                            star_mass_refined = rp_stars_halfmass[('stars', 'particle_mass')][halfmass_refined_bins] + stars_halfmass_to_add

                            stellar_halfmass_010_difference_refined.append(stellar_mass_010/2 - star_mass_refined)
                        print(stellar_halfmass_010_difference_refined)
                        for random_index in range(len(stellar_halfmass_010_difference_refined)):
                            print(random_index, stellar_halfmass_010_difference_refined[random_index])
                            if random_index == range(len(stellar_halfmass_010_difference_refined))[-1]:
                                print('Did not find stellar halfmass 010')
                                break
                            else:
                                if stellar_halfmass_010_difference_refined[random_index]*stellar_halfmass_010_difference_refined[random_index+1] < 0:
                                    radius_x = rp_stars_halfmass.x_bins[random_index]
                                    radius_y = rp_stars_halfmass.x_bins[random_index + 1]
                                    stellar_halfmass_010 = (radius_x + radius_y) / 2
                                    stellar_halfmass_010_delta = stellar_halfmass_010 - radius_x
                                    print('fixed 010')
                                    print(radius_x, rp_stars_halfmass[('stars', 'particle_mass')][random_index-1] + stars_halfmass_to_add)
                                    print(radius_y, rp_stars_halfmass[('stars', 'particle_mass')][random_index] + stars_halfmass_to_add)
                                    print(stellar_mass_010.in_units('Msun')/2)
                                    break

                    
                    if stellar_halfmass_015_difference[index_4] * stellar_halfmass_015_difference[index_4+1] < 0:
                        print('Finding 015 halfmass')
                        if index_4 == 0:
                            stars_halfmass_to_add = 0
                            print('0')
                        else:
                            stars_halfmass_to_add = rp_stars[('stars', 'particle_mass')][index_4-1]
                            print(stars_halfmass_to_add.in_units('g'))
                        
                        radius_3 = rp_stars.x_bins[index_4+0]
                        radius_4 = rp_stars.x_bins[index_4+1]
                        
                        radius_bins_refined_halfmass = ds.arr(np.linspace(radius_3, radius_4, num = 100), 'kpc')

                        rp_stars_halfmass = yt.create_profile(sp, ('stars', 'particle_radius'), [('stars', 'particle_mass')], accumulation = True,
                                                                 units = {('stars', 'particle_radius'): 'kpc', ('stars', 'particle_mass'): 'Msun'},
                                                                 weight_field=None, override_bins={('stars', 'particle_radius'): radius_bins_refined_halfmass})

                        stellar_halfmass_015_difference_refined = []

                        for halfmass_refined_bins in range(len(rp_stars_halfmass[('stars', 'particle_mass')])):
                            #star_radius_refined = rp_stars_halfmass.x_bins[halfmass_refined_bins + 1]
                            if halfmass_refined_bins == 0:
                                stellar_halfmass_015_difference_refined.append(stellar_mass_015/2 - stars_halfmass_to_add)
                                
                            star_mass_refined = rp_stars_halfmass[('stars', 'particle_mass')][halfmass_refined_bins] + stars_halfmass_to_add

                            stellar_halfmass_015_difference_refined.append(stellar_mass_015/2 - star_mass_refined)
                        print(stellar_halfmass_015_difference_refined)
                        for random_index in range(len(stellar_halfmass_015_difference_refined)):
                            print(random_index, stellar_halfmass_015_difference_refined[random_index])
                            if random_index == range(len(stellar_halfmass_015_difference_refined))[-1]:
                                print('Did not find stellar halfmass 015')
                                break
                            else:
                                if stellar_halfmass_015_difference_refined[random_index]*stellar_halfmass_015_difference_refined[random_index+1] < 0:
                                    radius_x = rp_stars_halfmass.x_bins[random_index]
                                    radius_y = rp_stars_halfmass.x_bins[random_index + 1]
                                    stellar_halfmass_015 = (radius_x + radius_y) / 2
                                    stellar_halfmass_015_delta = stellar_halfmass_015 - radius_x
                                    print('fixed 015')
                                    print(radius_x, rp_stars_halfmass[('stars', 'particle_mass')][random_index-1] + stars_halfmass_to_add)
                                    print(radius_y, rp_stars_halfmass[('stars', 'particle_mass')][random_index] + stars_halfmass_to_add)
                                    print(stellar_mass_015.in_units('Msun')/2)
                                    break

                    if stellar_halfmass_020_difference[index_4] * stellar_halfmass_020_difference[index_4+1] < 0:
                        print('Finding 020 halfmass')
                        if index_4 == 0:
                            stars_halfmass_to_add = 0
                            print('0')
                        else:
                            stars_halfmass_to_add = rp_stars[('stars', 'particle_mass')][index_4-1]
                            print(stars_halfmass_to_add.in_units('g'))
                        
                        radius_3 = rp_stars.x_bins[index_4+0]
                        radius_4 = rp_stars.x_bins[index_4+1]
                        
                        radius_bins_refined_halfmass = ds.arr(np.linspace(radius_3, radius_4, num = 100), 'kpc')

                        rp_stars_halfmass = yt.create_profile(sp, ('stars', 'particle_radius'), [('stars', 'particle_mass')], accumulation = True,
                                                                 units = {('stars', 'particle_radius'): 'kpc', ('stars', 'particle_mass'): 'Msun'},
                                                                 weight_field=None, override_bins={('stars', 'particle_radius'): radius_bins_refined_halfmass})

                        stellar_halfmass_020_difference_refined = []

                        for halfmass_refined_bins in range(len(rp_stars_halfmass[('stars', 'particle_mass')])):
                            #star_radius_refined = rp_stars_halfmass.x_bins[halfmass_refined_bins + 1]
                            if halfmass_refined_bins == 0:
                                stellar_halfmass_020_difference_refined.append(stellar_mass_020/2 - stars_halfmass_to_add)
                                
                            star_mass_refined = rp_stars_halfmass[('stars', 'particle_mass')][halfmass_refined_bins] + stars_halfmass_to_add

                            stellar_halfmass_020_difference_refined.append(stellar_mass_020/2 - star_mass_refined)
                        print(stellar_halfmass_020_difference_refined)
                        for random_index in range(len(stellar_halfmass_020_difference_refined)):
                            print(random_index, stellar_halfmass_020_difference_refined[random_index])
                            if random_index == range(len(stellar_halfmass_020_difference_refined))[-1]:
                                print('Did not find stellar halfmass 020')
                                break
                            else:
                                if stellar_halfmass_020_difference_refined[random_index]*stellar_halfmass_020_difference_refined[random_index+1] < 0:
                                    radius_x = rp_stars_halfmass.x_bins[random_index]
                                    radius_y = rp_stars_halfmass.x_bins[random_index + 1]
                                    stellar_halfmass_020 = (radius_x + radius_y) / 2
                                    stellar_halfmass_020_delta = stellar_halfmass_020 - radius_x
                                    print('fixed 020')
                                    print(radius_x, rp_stars_halfmass[('stars', 'particle_mass')][random_index-1] + stars_halfmass_to_add)
                                    print(radius_y, rp_stars_halfmass[('stars', 'particle_mass')][random_index] + stars_halfmass_to_add)
                                    print(stellar_mass_020.in_units('Msun')/2)
                                    break


                    

            
            #now add the values calculated to the lists for saving to a file
            rockstar_rvir_list.append(float(rvir.in_units('kpc')))
            rockstar_mvir_list.append(float(mvir.in_units('Msun')))
            rockstar_id_list.append(rockstar_id)
        
            new_rvir_list.append(float(new_rvir.in_units('kpc')))
            new_rvir_delta_list.append(float(new_rvir_delta.in_units('kpc')))
        
            stellar_mass_rvir_list.append(float(stellar_mass_rvir.in_units('Msun')))
            gas_mass_rvir_list.append(float(gas_mass_rvir.in_units('Msun')))
            darkmatter_mass_rvir_list.append(float(darkmatter_mass_rvir.in_units('Msun')))
                                                   
            stellar_mass_010_list.append(float(stellar_mass_010.in_units('Msun')))
            gas_mass_010_list.append(float(gas_mass_010.in_units('Msun')))
            stellar_mass_015_list.append(float(stellar_mass_015.in_units('Msun')))
            gas_mass_015_list.append(float(stellar_mass_015.in_units('Msun')))
            stellar_mass_020_list.append(float(stellar_mass_020.in_units('Msun')))
            gas_mass_020_list.append(float(gas_mass_020.in_units('Msun')))
        
            stellar_halfmass_010_list.append(float(stellar_halfmass_010.in_units('kpc'))) 
            stellar_halfmass_010_delta_list.append(float(stellar_halfmass_010_delta.in_units('kpc')))
            stellar_halfmass_015_list.append(float(stellar_halfmass_015.in_units('kpc'))) 
            stellar_halfmass_015_delta_list.append(float(stellar_halfmass_015_delta.in_units('kpc')))
            stellar_halfmass_020_list.append(float(stellar_halfmass_020.in_units('kpc'))) 
            stellar_halfmass_020_delta_list.append(float(stellar_halfmass_020_delta.in_units('kpc')))
            
        #now write the values to a file
        
        ascii_data = [rockstar_id_list, rockstar_mvir_list, rockstar_rvir_list, new_rvir_list, new_rvir_delta_list,
                     stellar_mass_rvir_list, gas_mass_rvir_list, darkmatter_mass_rvir_list, stellar_mass_010_list,
                     gas_mass_010_list, stellar_mass_015_list, gas_mass_015_list, stellar_mass_020_list,
                     gas_mass_020_list, stellar_halfmass_010_list, stellar_halfmass_010_delta_list,
                     stellar_halfmass_015_list, stellar_halfmass_015_delta_list,
                     stellar_halfmass_020_list, stellar_halfmass_020_delta_list]
        ascii_names = ['Rockstar_Id[0]', 'Rockstar_Mvir[1]', 'Rockstar_Rvir[2]', 'New_Rvir_Calculation[3]', 'New_Rvir_error[4]',
                      'Stellar_Mass_1Rvir[5]', 'Gas_Mass_1Rvir[6]', 'Darkmatter_Mass_1Rvir[7]', 'Stellar_Mass_.1Rvir[8]',
                      'Gas_Mass_.1Rvir[9]', 'Stellar_Mass_.15Rvir[10]', 'Gas_Mass_.15Rvir[11]', 'Stellar_Mass_.2Rvir[12]',
                      'Gas_Mass_.2Rvir[13]', 'Stellar_Halfmass_Radius_.1Rvir[14]', 'Stellar_Halfmass_Radius_.1Rvir_error[15]',
                      'Stellar_Halfmass_Radius_.15Rvir[16]', 'Stellar_Halfmass_Radius_.15Rvir_error[17]',
                      'Stellar_Halfmass_Radius_.2Rvir[18]', 'Stellar_Halfmass_Radius_.2Rvir_error[19]']
        
        file_name = '%s/improvedrvir%s.ascii' % (out_dir, current_scale)
        data = Table(ascii_data, names=ascii_names)
        ascii.write(data, output=file_name, overwrite=True) 

        profile_data = [rockstar_id_list, rockstar_mvir_list, rockstar_rvir_list, new_rvir_list, new_rvir_delta_list,
                       profile_radii_list, gas_profile_list, star_profile_list, darkmatter_profile_list]

        #profile_names = ['Rockstar_Id', 'Rockstar_Mvir', 'Rockstar_Rvir', 'New_Rvir_Calculation', 'New_Rvir_error',
        #                'Profile_Radii', 'Gas_Profile', 'Star_Profile', 'Darkmatter_Profile']
        #profile_data = Table(profile_data, names=profile_names)
        #profile_file_name = '%s/improvedrvirprofiles%s.ascii' % (out_dir, current_scale)
        #ascii.write(profile_data, output=profile_file_name, overwrite=True)

