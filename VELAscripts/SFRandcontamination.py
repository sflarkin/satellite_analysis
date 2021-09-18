import yt
import numpy as np
from yt import YTQuantity

import glob
import os, argparse

from satellite_analysis.catalogreaders import consistentcatalogreader as consistent

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('propermvircalc_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--remove_subhalos', nargs='?', default='True')
    parser.add_argument('--halo_mass', nargs='?', default=10e+09)
    args = vars(parser.parse_args())
    return args


args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
propermvircalc_dir = args['propermvircalc_dir']
out_dir = args['out_dir']
remove_subhalos = 'True'
halo_mass = 10e+09

#Check if the output directory exists and if it does not, create it.
if not os.path.exists(out_dir):
    print('Creating output directory')
    os.makedirs(out_dir)

VELA_snaps = glob.glob(VELA_dir + '/10MpcBox*')
VELA_snaps.sort()
print(VELA_snaps)

completed = glob.glob('{}/SFRandcontamination*.ascii'.format(out_dir))

consistent.consistent_catalog_reader(input_dir, remove_subhalos='True', halo_mass=10e9)

for snap in reversed(VELA_snaps):
    ds = yt.load(snap)
    current_scale = str(ds.scale_factor)[2:5]
    
    file_name = '%s/SFRandcontamination%s.ascii' % (out_dir, current_scale)
    
    if file_name in completed:
        print('Already completed Snapshot')
        continue
    
    #first we find the matching file in the consistent trees catalog so we can extract the center of the halo
        
    print('Trying to find matching Consistent File for scale: {}'.format(current_scale))
    
    position = [pos for pos, loc in enumerate(consistent.consistent_file_index) if loc == current_scale]

    if position == [] or len(position) > 1:
        print('Could not find corresponding consistent trees data for file {}'.format(snap))
        continue
    else:
        
        #now if we find the right consistent trees file, we look for the proper mvir calc file, to find the
        #virial radius to find the values within
        proper_rvir_list = []
        proper_mvir_list = []
        proper_id_list = []
        stellar010_rvir_list = []
        stellar015_rvir_list = []
        stellar020_rvir_list = []
        propermvircalc_file = '{}/improvedrvir{}.ascii'.format(propermvircalc_dir, current_scale) 
        
        try:
            readfile = open(propermvircalc_file)
            lines = readfile.readlines()
            
            for line in lines:
                #ignore the first line, as it is only text explaing what each column is
                if line == lines[0]:
                    continue
                
                current_line = line.split()
                halo_id = current_line[0]
                improved_rvir = current_line[3]
                improved_mvir = float(current_line[5]) + float(current_line[6]) + float(current_line[7])
                proper_id_list.append(halo_id)
                proper_mvir_list.append(float(improved_mvir))
                proper_rvir_list.append(float(improved_rvir))
                stellar010_rvir_list.append(float(current_line[8]))
                stellar015_rvir_list.append(float(current_line[10]))
                stellar020_rvir_list.append(float(current_line[12]))
                
            
        except FileNotFoundError:
            print('Could not find ProperMvir Data for file {}'.format(propermvircalc_file))
            continue  
        
        print('Matching Scales Test:{} {}'.format(current_scale, consistent.consistent_file_index[position[0]]))
        
        #make the lists to hold the data to be written to files later
        
        count = 0
        
        
        halo_data = consistent.halo_data_sorted[position[0]]
        print('Checking to see if length of rvir is equal to number of halos')
        print(len(halo_data) == len(proper_rvir_list))

        
        #make the lists to save the info for
        halo_id_list = []
        rvir_list = []
        mvir_list = []
        
        stars0_SFR_010rvir_list, stars0_SFR_015rvir_list, stars0_SFR_020rvir_list = [], [], []
        stars1_SFR_010rvir_list, stars1_SFR_015rvir_list, stars1_SFR_020rvir_list = [], [], []
        stars2_SFR_010rvir_list, stars2_SFR_015rvir_list, stars2_SFR_020rvir_list = [], [], []
        stars3_SFR_010rvir_list, stars3_SFR_015rvir_list, stars3_SFR_020rvir_list = [], [], []
        stars4_SFR_010rvir_list, stars4_SFR_015rvir_list, stars4_SFR_020rvir_list = [], [], []
        
        stars0_SSFR_010rvir_list, stars0_SSFR_015rvir_list, stars0_SSFR_020rvir_list = [], [], []
        stars1_SSFR_010rvir_list, stars1_SSFR_015rvir_list, stars1_SSFR_020rvir_list = [], [], []
        stars2_SSFR_010rvir_list, stars2_SSFR_015rvir_list, stars2_SSFR_020rvir_list = [], [], []
        stars3_SSFR_010rvir_list, stars3_SSFR_015rvir_list, stars3_SSFR_020rvir_list = [], [], []
        stars4_SSFR_010rvir_list, stars4_SSFR_015rvir_list, stars4_SSFR_020rvir_list = [], [], []
        
        rvir1_mass_contamination_list = []
        rvir_number_contamination_list = []
        rvir1_mass_contamination_fraction_list = []
        rvir1_number_contamination_fraction_list = []
        rvir1_larger_darkmatter_present_list = []
        rvir010_number_contamination_list = []
        
        stars100Myr_SFR_010rvir_list = [] 
        stars100Myr_SFR_015rvir_list = [] 
        stars100Myr_SFR_020rvir_list = []  
            
        stars100Myr_SSFR_010rvir_list = [] 
        stars100Myr_SSFR_015rvir_list = []
        stars100Myr_SSFR_020rvir_list = []
        
        stars60Myr_SFR_010rvir_list = [] 
        stars60Myr_SFR_015rvir_list = []
        stars60Myr_SFR_020rvir_list = []

        stars60Myr_SSFR_010rvir_list = []
        stars60Myr_SSFR_015rvir_list = []
        stars60Myr_SSFR_020rvir_list = []

        #the tacchella numbers 

        Tacchella5Myr_SFR_010_final_average_list = []
        Tacchella5Myr_SFR_010_final_median_list = []
        Tacchella5Myr_SSFR_010_final_average_list = []
        Tacchella5Myr_SSFR_010_final_median_list = []

        Tacchella5Myr_SFR_015_final_average_list = []
        Tacchella5Myr_SFR_015_final_median_list = []
        Tacchella5Myr_SSFR_015_final_average_list = []
        Tacchella5Myr_SSFR_015_final_median_list = []

        Tacchella5Myr_SFR_020_final_average_list = []
        Tacchella5Myr_SFR_020_final_median_list = []
        Tacchella5Myr_SSFR_020_final_average_list = []
        Tacchella5Myr_SSFR_020_final_median_list = []

        Tacchella60Myr_SFR_010_final_average_list = []
        Tacchella60Myr_SFR_010_final_median_list = []
        Tacchella60Myr_SSFR_010_final_average_list = []
        Tacchella60Myr_SSFR_010_final_median_list = []

        Tacchella60Myr_SFR_015_final_average_list = []
        Tacchella60Myr_SFR_015_final_median_list = []
        Tacchella60Myr_SSFR_015_final_average_list = []
        Tacchella60Myr_SSFR_015_final_median_list = []

        Tacchella60Myr_SFR_020_final_average_list = []
        Tacchella60Myr_SFR_020_final_median_list = []
        Tacchella60Myr_SSFR_020_final_average_list = []
        Tacchella60Myr_SSFR_020_final_median_list = []

        Tacchella100Myr_SFR_010_final_average_list = []
        Tacchella100Myr_SFR_010_final_median_list = []
        Tacchella100Myr_SSFR_010_final_average_list = []
        Tacchella100Myr_SSFR_010_final_median_list = []

        Tacchella100Myr_SFR_015_final_average_list = []
        Tacchella100Myr_SFR_015_final_median_list = []
        Tacchella100Myr_SSFR_015_final_average_list = []
        Tacchella100Myr_SSFR_015_final_median_list = []

        Tacchella100Myr_SFR_020_final_average_list = []
        Tacchella100Myr_SFR_020_final_median_list = []
        Tacchella100Myr_SSFR_020_final_average_list = []
        Tacchella100Myr_SSFR_020_final_median_list = []

        
        
        for index in range(len(halo_data)):
            
            halo = halo_data[index]
            
            
            #the below if else doesnt work due to length of halo_data being longer in the case one was skipped
            #just do else always as it works for every case
            right_index = [pos for pos, loc in enumerate(proper_id_list) if loc == str(float(halo[1]))]
            if right_index == []:
            print('halo not found in propermvircalc')
                continue
            x = ds.quan(float(halo[17]), 'Mpccm/h')
            y = ds.quan(float(halo[18]), 'Mpccm/h')
            z = ds.quan(float(halo[19]), 'Mpccm/h')
            center = [x,y,z]
            rvir = proper_rvir_list[right_index[0]]
            mvir = proper_mvir_list[right_index[0]]
            Mstar010 = stellar010_rvir_list[right_index[0]]
            Mstar015 = stellar015_rvir_list[right_index[0]]
            Mstar020 = stellar020_rvir_list[right_index[0]]
            
            #print(halo[1], proper_id_list[index])
            #if float(halo[1]) == float(proper_id_list[index]):
                #x = ds.quan(float(halo[17]), 'Mpccm/h')
                #y = ds.quan(float(halo[18]), 'Mpccm/h')
                #z = ds.quan(float(halo[19]), 'Mpccm/h')
                #center = [x,y,z]
                #rvir = proper_rvir_list[index]
                #mvir = proper_mvir_list[index]
                #Mstar010 = stellar010_rvir_list[index]
                #Mstar015 = stellar015_rvir_list[index]
                #Mstar020 = stellar020_rvir_list[index]
        
            #so this doesnt work if we skipped the halo in propermvircalc
            #so we want to search to find the right index for proper lists
            #but only if they dont match normally
            #testing in FFW fits for profiles
            
            #else:
                #we try to enumerate to find the matching halo_id, if not found skip
                #right_index = [pos for pos, loc in enumerate(proper_id_list) if loc == str(float(halo[1]))]
                #if right_index == []:
                    #print('halo not found in propermvircalc')
                    #continue
                #x = ds.quan(float(halo[17]), 'Mpccm/h')
                #y = ds.quan(float(halo[18]), 'Mpccm/h')
                #z = ds.quan(float(halo[19]), 'Mpccm/h')
                #center = [x,y,z]
                #rvir = proper_rvir_list[right_index[0]]
                #mvir = proper_mvir_list[right_index[0]]
                #Mstar010 = stellar010_rvir_list[right_index[0]]
                #Mstar015 = stellar015_rvir_list[right_index[0]]
                #Mstar020 = stellar020_rvir_list[right_index[0]]
            
            ad = ds.all_data()
            
            rvir1_sphere = ds.sphere(center, (rvir, 'kpc'))
            rvir010_sphere = ds.sphere(center, (0.1*rvir, 'kpc'))
            rvir015_sphere = ds.sphere(center, (0.15*rvir, 'kpc'))
            rvir020_sphere = ds.sphere(center, (0.2*rvir, 'kpc'))
            
            
            #now set up the filters to find the data we want in the regions
            
            #we need a filter for the 5 stellar populations by age
            #use each of these on the .1 .15 and .2 rvir spheres, so 15 total spheres
            star_creation_times = yt.np.unique(ad[('stars', 'particle_creation_time')])
            def youngest_one_stars(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[-1]
                return filter
            yt.add_particle_filter('stars0', function=youngest_one_stars, filtered_type='stars', requires=['particle_creation_time'])
            ds.add_particle_filter('stars0')
           
            print('success') 
            def youngest_two_stars(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[-2]
                return filter
            yt.add_particle_filter('stars1', function=youngest_two_stars, filtered_type='stars', requires=['particle_creation_time'])
            ds.add_particle_filter('stars1')
            
            def youngest_three_stars(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[-3]
                return filter
            yt.add_particle_filter('stars2', function=youngest_three_stars, filtered_type='stars', requires=['particle_creation_time'])
            ds.add_particle_filter('stars2')
            
            def youngest_four_stars(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[-4]
                return filter
            yt.add_particle_filter('stars3', function=youngest_four_stars, filtered_type='stars', requires=['particle_creation_time'])
            ds.add_particle_filter('stars3')
            
            def youngest_five_stars(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[-5]
                return filter
            yt.add_particle_filter('stars4', function=youngest_five_stars, filtered_type='stars', requires=['particle_creation_time'])
            ds.add_particle_filter('stars4')
            
            
            current_time = ds.current_time.in_units('yr')
            print(current_time)
            stars0_time = yt.np.unique(ad[('stars0', 'particle_creation_time')])[0].in_units('yr')
            stars1_time = yt.np.unique(ad[('stars1', 'particle_creation_time')])[0].in_units('yr')
            stars2_time = yt.np.unique(ad[('stars2', 'particle_creation_time')])[0].in_units('yr')
            stars3_time = yt.np.unique(ad[('stars3', 'particle_creation_time')])[0].in_units('yr')
            stars4_time = yt.np.unique(ad[('stars4', 'particle_creation_time')])[0].in_units('yr')
            
            #find the mass of the population for each filter for each sphere
            
    
            #first for the .1 Rvir
            stars0_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars0', 'particle_mass_initial')]).in_units('Msun')
            stars1_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars1', 'particle_mass_initial')]).in_units('Msun')
            stars2_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars2', 'particle_mass_initial')]).in_units('Msun')
            stars3_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars3', 'particle_mass_initial')]).in_units('Msun')
            stars4_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars4', 'particle_mass_initial')]).in_units('Msun')
            
            stars0_SFR_010rvir = float(stars0_mass_010rvir / (current_time - stars0_time))
            stars1_SFR_010rvir = float(stars1_mass_010rvir / (current_time - stars1_time))
            stars2_SFR_010rvir = float(stars2_mass_010rvir / (current_time - stars2_time))
            stars3_SFR_010rvir = float(stars3_mass_010rvir / (current_time - stars3_time))
            stars4_SFR_010rvir = float(stars4_mass_010rvir / (current_time - stars4_time))
            
            #get the Mstar for the SSFR
            #changed to get these earlier for propermvircalc skip fix
            
            stars0_SSFR_010rvir = stars0_SFR_010rvir/Mstar010
            stars1_SSFR_010rvir = stars1_SFR_010rvir/Mstar010
            stars2_SSFR_010rvir = stars2_SFR_010rvir/Mstar010
            stars3_SSFR_010rvir = stars3_SFR_010rvir/Mstar010
            stars4_SSFR_010rvir = stars4_SFR_010rvir/Mstar010
            
            print(stars0_SFR_010rvir)
            print(stars1_SFR_010rvir)
            print(stars2_SFR_010rvir)
            print(stars3_SFR_010rvir)
            print(stars4_SFR_010rvir)
            
            
            
            
            
            #now for .15rvir
            
            stars0_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars0', 'particle_mass_initial')]).in_units('Msun')
            stars1_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars1', 'particle_mass_initial')]).in_units('Msun')
            stars2_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars2', 'particle_mass_initial')]).in_units('Msun')
            stars3_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars3', 'particle_mass_initial')]).in_units('Msun')
            stars4_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars4', 'particle_mass_initial')]).in_units('Msun')
            
            stars0_SFR_015rvir = float(stars0_mass_015rvir / (current_time - stars0_time))
            stars1_SFR_015rvir = float(stars1_mass_015rvir / (current_time - stars1_time))
            stars2_SFR_015rvir = float(stars2_mass_015rvir / (current_time - stars2_time))
            stars3_SFR_015rvir = float(stars3_mass_015rvir / (current_time - stars3_time))
            stars4_SFR_015rvir = float(stars4_mass_015rvir / (current_time - stars4_time))
            
            stars0_SSFR_015rvir = stars0_SFR_015rvir/Mstar015
            stars1_SSFR_015rvir = stars1_SFR_015rvir/Mstar015
            stars2_SSFR_015rvir = stars2_SFR_015rvir/Mstar015
            stars3_SSFR_015rvir = stars3_SFR_015rvir/Mstar015
            stars4_SSFR_015rvir = stars4_SFR_015rvir/Mstar015

            
            
            
            #now for .2Rvir
            stars0_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars0', 'particle_mass_initial')]).in_units('Msun')
            stars1_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars1', 'particle_mass_initial')]).in_units('Msun')
            stars2_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars2', 'particle_mass_initial')]).in_units('Msun')
            stars3_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars3', 'particle_mass_initial')]).in_units('Msun')
            stars4_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars4', 'particle_mass_initial')]).in_units('Msun')
            
            stars0_SFR_020rvir = float(stars0_mass_020rvir / (current_time - stars0_time))
            stars1_SFR_020rvir = float(stars1_mass_020rvir / (current_time - stars1_time))
            stars2_SFR_020rvir = float(stars2_mass_020rvir / (current_time - stars2_time))
            stars3_SFR_020rvir = float(stars3_mass_020rvir / (current_time - stars3_time))
            stars4_SFR_020rvir = float(stars4_mass_020rvir / (current_time - stars4_time))
            
            stars0_SSFR_020rvir = stars0_SFR_020rvir/Mstar020
            stars1_SSFR_020rvir = stars1_SFR_020rvir/Mstar020
            stars2_SSFR_020rvir = stars2_SFR_020rvir/Mstar020
            stars3_SSFR_020rvir = stars3_SFR_020rvir/Mstar020
            stars4_SSFR_020rvir = stars4_SFR_020rvir/Mstar020
            
            
            #now we do the 100Myr one if we can

            count = -1
            Myr100_time_count = None
            Myr100 = ds.quan(100, 'Myr')
            while count >= -len(star_creation_times):
                stars100Myr_time_candidate = (current_time - star_creation_times[count])
                if stars100Myr_time_candidate >= Myr100:
                    Myr100_time_count = count
                    break
                count = count - 1

            if Myr100_time_count == None:
                print('Not enough timesteps for 100 Myr to pass')

                #set the values to 0
                stars100Myr_SFR_010rvir = 0 
                stars100Myr_SFR_015rvir = 0 
                stars100Myr_SFR_020rvir = 0 

                #SSFR
                stars100Myr_SSFR_010rvir = 0
                stars100Myr_SSFR_015rvir = 0
                stars100Myr_SSFR_020rvir = 0

            else:
                def Myr100_stars(pfilter, data):
                    filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[Myr100_time_count]
                    return filter
                yt.add_particle_filter('stars100Myr', function=Myr100_stars, filtered_type='stars', requires=['particle_creation_time'])
                ds.add_particle_filter('stars100Myr')


                stars100Myr_time = yt.np.unique(ad[('stars100Myr', 'particle_creation_time')])[0].in_units('yr')

                #SFR
                stars100Myr_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars100Myr', 'particle_mass_initial')]).in_units('Msun')
                stars100Myr_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars100Myr', 'particle_mass_initial')]).in_units('Msun')
                stars100Myr_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars100Myr', 'particle_mass_initial')]).in_units('Msun')
                
                stars100Myr_SFR_010rvir = float(stars100Myr_mass_010rvir / (current_time - stars100Myr_time).in_units('yr')) 
                stars100Myr_SFR_015rvir = float(stars100Myr_mass_015rvir / (current_time - stars100Myr_time).in_units('yr')) 
                stars100Myr_SFR_020rvir = float(stars100Myr_mass_020rvir / (current_time - stars100Myr_time).in_units('yr')) 

                #SSFR
                stars100Myr_SSFR_010rvir = stars100Myr_SFR_010rvir/Mstar010
                stars100Myr_SSFR_015rvir = stars100Myr_SFR_015rvir/Mstar015
                stars100Myr_SSFR_020rvir = stars100Myr_SFR_020rvir/Mstar020

            
            
            
            #now we do the 60Myr one if we can
            

            count = -1
            Myr60_time_count = None
            Myr60 = ds.quan(60, 'Myr')
            while count >= -len(star_creation_times):
                stars60Myr_time_candidate = (current_time - star_creation_times[count])
                if stars60Myr_time_candidate >= Myr60:
                    Myr60_time_count = count
                    break
                count = count - 1

            if Myr60_time_count == None:
                print('Not enough timesteps for 100 Myr to pass')

                #set the values to 0
                stars60Myr_SFR_010rvir = 0 
                stars60Myr_SFR_015rvir = 0 
                stars60Myr_SFR_020rvir = 0 

                #SSFR
                stars60Myr_SSFR_010rvir = 0
                stars60Myr_SSFR_015rvir = 0
                stars60Myr_SSFR_020rvir = 0

            else:
                def Myr60_stars(pfilter, data):
                    filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= star_creation_times[Myr60_time_count]
                    return filter
                yt.add_particle_filter('stars60Myr', function=Myr60_stars, filtered_type='stars', requires=['particle_creation_time'])
                ds.add_particle_filter('stars60Myr')


                stars60Myr_time = yt.np.unique(ad[('stars60Myr', 'particle_creation_time')])[0].in_units('yr')

                #SFR
                stars60Myr_mass_010rvir = rvir010_sphere.quantities.total_quantity([('stars60Myr', 'particle_mass_initial')]).in_units('Msun')
                stars60Myr_mass_015rvir = rvir015_sphere.quantities.total_quantity([('stars60Myr', 'particle_mass_initial')]).in_units('Msun')
                stars60Myr_mass_020rvir = rvir020_sphere.quantities.total_quantity([('stars60Myr', 'particle_mass_initial')]).in_units('Msun')
                
                stars60Myr_SFR_010rvir = float(stars60Myr_mass_010rvir / (current_time - stars60Myr_time).in_units('yr')) 
                stars60Myr_SFR_015rvir = float(stars60Myr_mass_015rvir / (current_time - stars60Myr_time).in_units('yr')) 
                stars60Myr_SFR_020rvir = float(stars60Myr_mass_020rvir / (current_time - stars60Myr_time).in_units('yr')) 

                #SSFR
                stars60Myr_SSFR_010rvir = stars60Myr_SFR_010rvir/Mstar010
                stars60Myr_SSFR_015rvir = stars60Myr_SFR_015rvir/Mstar015
                stars60Myr_SSFR_020rvir = stars60Myr_SFR_020rvir/Mstar020
            
            
            #these are the tacchella methods of SFR for his paper (60 Myr), as well as 5 and 100 Myr for comparison 

            Myr10 = ds.quan(10, 'Myr')
            Myr80 = ds.quan(80, 'Myr')
            Myr120 = ds.quan(120, 'Myr')

            if (ds.current_time - star_creation_times[0]) > Myr80:
                #this means there are enough timesteps to do this
                age_cuts_60Myr = np.arange(40, 80.2, 0.2)

                Tacchella60Myr_SFR_010rvir_list = []
                Tacchella60Myr_SFR_015rvir_list = []
                Tacchella60Myr_SFR_020rvir_list = []

                for age in age_cuts_60Myr:

                    current_age = ds.quan(age, 'Myr')

                    def Tacchella_60Myr_stars(pfilter, data):
                        filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= (ds.current_time - current_age)
                        return filter
                    yt.add_particle_filter('Tacchellastars60Myr', function=Tacchella_60Myr_stars, filtered_type='stars', requires=['particle_creation_time'])
                    ds.add_particle_filter('Tacchellastars60Myr')


                    Tacchella60Myr_mass_010rvir = rvir010_sphere.quantities.total_quantity([('Tacchellastars60Myr', 'particle_mass_initial')]).in_units('Msun')
                    Tacchella60Myr_mass_015rvir = rvir015_sphere.quantities.total_quantity([('Tacchellastars60Myr', 'particle_mass_initial')]).in_units('Msun')
                    Tacchella60Myr_mass_020rvir = rvir020_sphere.quantities.total_quantity([('Tacchellastars60Myr', 'particle_mass_initial')]).in_units('Msun')

                    Tacchella60Myr_SFR_010rvir = float(Tacchella60Myr_mass_010rvir / (current_age.in_units('yr'))) 
                    Tacchella60Myr_SFR_015rvir = float(Tacchella60Myr_mass_015rvir / (current_age.in_units('yr'))) 
                    Tacchella60Myr_SFR_020rvir = float(Tacchella60Myr_mass_020rvir / (current_age.in_units('yr'))) 

                    Tacchella60Myr_SFR_010rvir_list.append(Tacchella60Myr_SFR_010rvir)
                    Tacchella60Myr_SFR_015rvir_list.append(Tacchella60Myr_SFR_015rvir)
                    Tacchella60Myr_SFR_020rvir_list.append(Tacchella60Myr_SFR_020rvir)


                #now we find the median and mean for the lists of all the values

                array_010 = np.asarray(Tacchella60Myr_SFR_010rvir_list)
                Tacchella60Myr_SFR_010_final_average = np.average(array_010)
                Tacchella60Myr_SFR_010_final_median = np.median(array_010)
                Tacchella60Myr_SSFR_010_final_average = Tacchella60Myr_SFR_010_final_average/Mstar010
                Tacchella60Myr_SSFR_010_final_median = Tacchella60Myr_SFR_010_final_median/Mstar010

                array_015 = np.asarray(Tacchella60Myr_SFR_015rvir_list)
                Tacchella60Myr_SFR_015_final_average = np.average(array_015)
                Tacchella60Myr_SFR_015_final_median = np.median(array_015)
                Tacchella60Myr_SSFR_015_final_average = Tacchella60Myr_SFR_015_final_average/Mstar015
                Tacchella60Myr_SSFR_015_final_median = Tacchella60Myr_SFR_015_final_median/Mstar015

                array_020 = np.asarray(Tacchella60Myr_SFR_020rvir_list)
                Tacchella60Myr_SFR_020_final_average = np.average(array_020)
                Tacchella60Myr_SFR_020_final_median = np.median(array_020)
                Tacchella60Myr_SSFR_020_final_average = Tacchella60Myr_SFR_020_final_average/Mstar020
                Tacchella60Myr_SSFR_020_final_median = Tacchella60Myr_SFR_020_final_median/Mstar020


            else:
                #set the variables for SFR and SSFR for these to 0

                Tacchella60Myr_SFR_010_final_average = 0
                Tacchella60Myr_SFR_010_final_median = 0
                Tacchella60Myr_SSFR_010_final_average = 0
                Tacchella60Myr_SSFR_010_final_median = 0

                Tacchella60Myr_SFR_015_final_average = 0
                Tacchella60Myr_SFR_015_final_median = 0
                Tacchella60Myr_SSFR_015_final_average = 0
                Tacchella60Myr_SSFR_015_final_median = 0

                Tacchella60Myr_SFR_020_final_average = 0
                Tacchella60Myr_SFR_020_final_median = 0
                Tacchella60Myr_SSFR_020_final_average = 0
                Tacchella60Myr_SSFR_020_final_median = 0



            if (ds.current_time - star_creation_times[0]) > Myr10:
                #this means there are enough timesteps to do this
                age_cuts_5Myr = np.arange(0.2, 10.2, 0.2)

                Tacchella5Myr_SFR_010rvir_list = []
                Tacchella5Myr_SFR_015rvir_list = []
                Tacchella5Myr_SFR_020rvir_list = []

                for age in age_cuts_5Myr:

                    current_age = ds.quan(age, 'Myr')

                    def Tacchella_5Myr_stars(pfilter, data):
                        filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= (ds.current_time - current_age)
                        return filter
                    yt.add_particle_filter('Tacchellastars5Myr', function=Tacchella_5Myr_stars, filtered_type='stars', requires=['particle_creation_time'])
                    ds.add_particle_filter('Tacchellastars5Myr')


                    Tacchella5Myr_mass_010rvir = rvir010_sphere.quantities.total_quantity([('Tacchellastars5Myr', 'particle_mass_initial')]).in_units('Msun')
                    Tacchella5Myr_mass_015rvir = rvir015_sphere.quantities.total_quantity([('Tacchellastars5Myr', 'particle_mass_initial')]).in_units('Msun')
                    Tacchella5Myr_mass_020rvir = rvir020_sphere.quantities.total_quantity([('Tacchellastars5Myr', 'particle_mass_initial')]).in_units('Msun')

                    Tacchella5Myr_SFR_010rvir = float(Tacchella5Myr_mass_010rvir / (current_age.in_units('yr'))) 
                    Tacchella5Myr_SFR_015rvir = float(Tacchella5Myr_mass_015rvir / (current_age.in_units('yr'))) 
                    Tacchella5Myr_SFR_020rvir = float(Tacchella5Myr_mass_020rvir / (current_age.in_units('yr'))) 

                    Tacchella5Myr_SFR_010rvir_list.append(Tacchella5Myr_SFR_010rvir)
                    Tacchella5Myr_SFR_015rvir_list.append(Tacchella5Myr_SFR_015rvir)
                    Tacchella5Myr_SFR_020rvir_list.append(Tacchella5Myr_SFR_020rvir)


                #now we find the median and mean for the lists of all the values
                
                array_010 = np.asarray(Tacchella5Myr_SFR_010rvir_list)
                Tacchella5Myr_SFR_010_final_average = np.average(array_010)
                Tacchella5Myr_SFR_010_final_median = np.median(array_010)
                Tacchella5Myr_SSFR_010_final_average = Tacchella5Myr_SFR_010_final_average/Mstar010
                Tacchella5Myr_SSFR_010_final_median = Tacchella5Myr_SFR_010_final_median/Mstar010

                array_015 = np.asarray(Tacchella5Myr_SFR_015rvir_list)
                Tacchella5Myr_SFR_015_final_average = np.average(array_015)
                Tacchella5Myr_SFR_015_final_median = np.median(array_015)
                Tacchella5Myr_SSFR_015_final_average = Tacchella5Myr_SFR_015_final_average/Mstar015
                Tacchella5Myr_SSFR_015_final_median = Tacchella5Myr_SFR_015_final_median/Mstar015

                array_020 = np.asarray(Tacchella5Myr_SFR_020rvir_list)
                Tacchella5Myr_SFR_020_final_average = np.average(array_020)
                Tacchella5Myr_SFR_020_final_median = np.median(array_020)
                Tacchella5Myr_SSFR_020_final_average = Tacchella5Myr_SFR_020_final_average/Mstar020
                Tacchella5Myr_SSFR_020_final_median = Tacchella5Myr_SFR_020_final_median/Mstar020
                
                print(array_010)
                print(array_015)
                print(array_020)

            else:
                #set the variables for SFR and SSFR for these to 0

                Tacchella5Myr_SFR_010_final_average = 0
                Tacchella5Myr_SFR_010_final_median = 0
                Tacchella5Myr_SSFR_010_final_average = 0
                Tacchella5Myr_SSFR_010_final_median = 0

                Tacchella5Myr_SFR_015_final_average = 0
                Tacchella5Myr_SFR_015_final_median = 0
                Tacchella5Myr_SSFR_015_final_average = 0
                Tacchella5Myr_SSFR_015_final_median = 0

                Tacchella5Myr_SFR_020_final_average = 0
                Tacchella5Myr_SFR_020_final_median = 0
                Tacchella5Myr_SSFR_020_final_average = 0
                Tacchella5Myr_SSFR_020_final_median = 0



            if (ds.current_time - star_creation_times[0]) > Myr120:
                #this means there are enough timesteps to do this
                age_cuts_100Myr = np.arange(80, 120.2, 0.2)

                Tacchella100Myr_SFR_010rvir_list = []
                Tacchella100Myr_SFR_015rvir_list = []
                Tacchella100Myr_SFR_020rvir_list = []

                for age in age_cuts_100Myr:

                    current_age = ds.quan(age, 'Myr')

                    def Tacchella_100Myr_stars(pfilter, data):
                        filter = data[(pfilter.filtered_type, 'particle_creation_time')] >= (ds.current_time - current_age)
                        return filter
                    yt.add_particle_filter('Tacchellastars100Myr', function=Tacchella_100Myr_stars, filtered_type='stars', requires=['particle_creation_time'])
                    ds.add_particle_filter('Tacchellastars100Myr')


                    Tacchella100Myr_mass_010rvir = rvir010_sphere.quantities.total_quantity([('Tacchellastars100Myr', 'particle_mass_initial')]).in_units('Msun')
                    Tacchella100Myr_mass_015rvir = rvir015_sphere.quantities.total_quantity([('Tacchellastars100Myr', 'particle_mass_initial')]).in_units('Msun')
                    Tacchella100Myr_mass_020rvir = rvir020_sphere.quantities.total_quantity([('Tacchellastars100Myr', 'particle_mass_initial')]).in_units('Msun')

                    Tacchella100Myr_SFR_010rvir = float(Tacchella100Myr_mass_010rvir / (current_age.in_units('yr'))) 
                    Tacchella100Myr_SFR_015rvir = float(Tacchella100Myr_mass_015rvir / (current_age.in_units('yr'))) 
                    Tacchella100Myr_SFR_020rvir = float(Tacchella100Myr_mass_020rvir / (current_age.in_units('yr'))) 

                    Tacchella100Myr_SFR_010rvir_list.append(Tacchella100Myr_SFR_010rvir)
                    Tacchella100Myr_SFR_015rvir_list.append(Tacchella100Myr_SFR_015rvir)
                    Tacchella100Myr_SFR_020rvir_list.append(Tacchella100Myr_SFR_020rvir)


                #now we find the median and mean for the lists of all the values

                array_010 = np.asarray(Tacchella100Myr_SFR_010rvir_list)
                Tacchella100Myr_SFR_010_final_average = np.average(array_010)
                Tacchella100Myr_SFR_010_final_median = np.median(array_010)
                Tacchella100Myr_SSFR_010_final_average = Tacchella100Myr_SFR_010_final_average/Mstar010
                Tacchella100Myr_SSFR_010_final_median = Tacchella100Myr_SFR_010_final_median/Mstar010

                array_015 = np.asarray(Tacchella100Myr_SFR_015rvir_list)
                Tacchella100Myr_SFR_015_final_average = np.average(array_015)
                Tacchella100Myr_SFR_015_final_median = np.median(array_015)
                Tacchella100Myr_SSFR_015_final_average = Tacchella100Myr_SFR_015_final_average/Mstar015
                Tacchella100Myr_SSFR_015_final_median = Tacchella100Myr_SFR_015_final_median/Mstar015

                array_020 = np.asarray(Tacchella100Myr_SFR_020rvir_list)
                Tacchella100Myr_SFR_020_final_average = np.average(array_020)
                Tacchella100Myr_SFR_020_final_median = np.median(array_020)
                Tacchella100Myr_SSFR_020_final_average = Tacchella100Myr_SFR_020_final_average/Mstar020
                Tacchella100Myr_SSFR_020_final_median = Tacchella100Myr_SFR_020_final_median/Mstar020


            else:
                #set the variables for SFR and SSFR for these to 0

                Tacchella100Myr_SFR_010_final_average = 0
                Tacchella100Myr_SFR_010_final_median = 0
                Tacchella100Myr_SSFR_010_final_average = 0
                Tacchella100Myr_SSFR_010_final_median = 0

                Tacchella100Myr_SFR_015_final_average = 0
                Tacchella100Myr_SFR_015_final_median = 0
                Tacchella100Myr_SSFR_015_final_average = 0
                Tacchella100Myr_SSFR_015_final_median = 0

                Tacchella100Myr_SFR_020_final_average = 0
                Tacchella100Myr_SFR_020_final_median = 0
                Tacchella100Myr_SSFR_020_final_average = 0
                Tacchella100Myr_SSFR_020_final_median = 0

            
            
            
            #we need a filter for the larger darkmatter particles, and for only the stage one up to see if they
            #are the same
            
            darkmatter_masses = yt.np.unique(ad[('darkmatter', 'particle_mass')])
            
            def darkmatter_mass_filter_all_larger(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_mass')] > darkmatter_masses[0]
                return filter
            yt.add_particle_filter('darkmatter0', function=darkmatter_mass_filter_all_larger, filtered_type='darkmatter', requires=['particle_mass'])
            ds.add_particle_filter('darkmatter0')
            
            def darkmatter_mass_filter_second_smallest(pfilter, data):
                filter = data[(pfilter.filtered_type, 'particle_mass')] == darkmatter_masses[1]
                return filter
            yt.add_particle_filter('darkmatter1', darkmatter_mass_filter_second_smallest, filtered_type='darkmatter', requires=['particle_mass'])
            ds.add_particle_filter('darkmatter1')
            
            
            rvir1_mass_contamination = float(rvir1_sphere.quantities.total_quantity([('darkmatter0', 'particle_mass')]).in_units('Msun'))
            rvir1_number_contamination = len(rvir1_sphere[('darkmatter0', 'particle_mass')])
            
            rvir1_mass_contamination_fraction = rvir1_mass_contamination/mvir
            rvir1_number_contamination_fraction = rvir1_number_contamination/len(rvir1_sphere[('darkmatter', 'particle_mass')])
            rvir1_larger_darkmatter_present = (rvir1_number_contamination == len(rvir1_sphere[('darkmatter1', 'particle_mass')]))
            
            #if contamination is above 1%, plot where the contamination is
            print(rvir1_mass_contamination_fraction)
            if rvir1_mass_contamination_fraction > 0.01:
                print('Contamination above 1%, plotting image')
                
                fig = plt.figure()
    
                grid = AxesGrid(fig, (0.075,0.075,15,5),
                nrows_ncols = (2, 3),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")
            
                a = yt.ParticlePlot(ds, ('stars', 'particle_position_x'), ('stars', 'particle_position_y'),\
                                      ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                b = yt.ParticlePlot(ds, ('stars', 'particle_position_y'), ('stars', 'particle_position_z'),\
                                      ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                c = yt.ParticlePlot(ds, ('stars', 'particle_position_z'), ('stars', 'particle_position_x'),\
                                      ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
            
                d = yt.ParticlePlot(ds, ('darkmatter0', 'particle_position_x'), ('darkmatter0', 'particle_position_y'),\
                                      ('darkmatter0', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                e = yt.ParticlePlot(ds, ('darkmatter0', 'particle_position_y'), ('darkmatter0', 'particle_position_z'),\
                                      ('darkmatter0', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                f = yt.ParticlePlot(ds, ('darkmatter0', 'particle_position_z'), ('darkmatter0', 'particle_position_x'),\
                                      ('darkmatter0', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
            
                a.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                b.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                c.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                d.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                e.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                f.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            
                index1 = 0
                for letter in [a,b,c]:
                    plot = letter.plots[('stars', 'particle_mass')]
                    plot.figure = fig
                    plot.axes = grid[index1].axes
                    plot.cax = grid.cbar_axes[0]
                    letter._setup_plots()
                    index1 = index1 + 1

                for letter in [d,e,f]:
                    plot = letter.plots[('darkmatter0', 'particle_mass')]
                    plot.figure = fig
                    plot.axes = grid[index1].axes
                    plot.cax = grid.cbar_axes[0]
                    letter._setup_plots()
                    index1 = index1 + 1
            
                plt.savefig('{}/darkmatter_contanimation_a{}_id{}.png'.format(out_dir, current_scale, halo[1]), bbox_inches='tight')
                plt.close()
            
            #if there are dark matter particles larger than the second smallest, plot where they are
            if rvir1_number_contamination != len(rvir1_sphere[('darkmatter1', 'particle_mass')]):
                print(rvir1_number_contamination, len(rvir1_sphere[('darkmatter1', 'particle_mass')]))
                print('Some Third smalles or larger particles present within Rvir')
                print('Generating Graph of Their locations')
                
                def darkmatter_mass_filter_larger_second_smallest(pfilter, data):
                    filter = data[(pfilter.filtered_type, 'particle_mass')] > darkmatter_masses[1]
                    return filter
                yt.add_particle_filter('darkmatter2', darkmatter_mass_filter_larger_second_smallest, filtered_type='darkmatter', requires=['particle_mass'])
                ds.add_particle_filter('darkmatter2')
                
                fig = plt.figure()
    
                grid = AxesGrid(fig, (0.075,0.075,15,5),
                nrows_ncols = (2, 3),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")
            
                a = yt.ParticlePlot(ds, ('stars', 'particle_position_x'), ('stars', 'particle_position_y'),\
                                      ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                b = yt.ParticlePlot(ds, ('stars', 'particle_position_y'), ('stars', 'particle_position_z'),\
                                      ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                c = yt.ParticlePlot(ds, ('stars', 'particle_position_z'), ('stars', 'particle_position_x'),\
                                      ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
            
                d = yt.ParticlePlot(ds, ('darkmatter2', 'particle_position_x'), ('darkmatter2', 'particle_position_y'),\
                                      ('darkmatter2', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                e = yt.ParticlePlot(ds, ('darkmatter2', 'particle_position_y'), ('darkmatter2', 'particle_position_z'),\
                                      ('darkmatter2', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
                f = yt.ParticlePlot(ds, ('darkmatter2', 'particle_position_z'), ('darkmatter2', 'particle_position_x'),\
                                      ('darkmatter2', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=rvir1_sphere)
            
                a.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                b.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                c.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                d.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                e.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                f.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            
                index1 = 0
                for letter in [a,b,c]:
                    plot = letter.plots[('stars', 'particle_mass')]
                    plot.figure = fig
                    plot.axes = grid[index1].axes
                    plot.cax = grid.cbar_axes[0]
                    letter._setup_plots()
                    index1 = index1 + 1

                for letter in [d,e,f]:
                    plot = letter.plots[('darkmatter2', 'particle_mass')]
                    plot.figure = fig
                    plot.axes = grid[index1].axes
                    plot.cax = grid.cbar_axes[0]
                    letter._setup_plots()
                    index1 = index1 + 1
            
                plt.savefig('{}/darkmatter_contanimation_larger_particles{}_a{}_id{}.png'.format(out_dir, rvir1_number_contamination-len(rvir1_sphere[('darkmatter1', 'particle_mass')]), current_scale, halo[1]), bbox_inches='tight')
                plt.close()
            
            #now see if there is contamination within .1Rvir
            
            rvir010_number_contamination = len(rvir010_sphere[('darkmatter0', 'particle_mass')])
            
            #now save the values to the lists for storing later
            
            
            
            halo_id_list.append(halo[1])
            rvir_list.append(rvir)
            mvir_list.append(mvir)
        
            stars0_SFR_010rvir_list.append(stars0_SFR_010rvir)
            stars0_SFR_015rvir_list.append(stars0_SFR_015rvir)
            stars0_SFR_020rvir_list.append(stars0_SFR_020rvir)
            stars1_SFR_010rvir_list.append(stars1_SFR_010rvir)
            stars1_SFR_015rvir_list.append(stars1_SFR_015rvir)
            stars1_SFR_020rvir_list.append(stars1_SFR_020rvir)
            stars2_SFR_010rvir_list.append(stars2_SFR_010rvir)
            stars2_SFR_015rvir_list.append(stars2_SFR_015rvir)
            stars2_SFR_020rvir_list.append(stars2_SFR_020rvir)
            stars3_SFR_010rvir_list.append(stars3_SFR_010rvir)
            stars3_SFR_015rvir_list.append(stars3_SFR_015rvir)
            stars3_SFR_020rvir_list.append(stars3_SFR_020rvir)
            stars4_SFR_010rvir_list.append(stars4_SFR_010rvir)
            stars4_SFR_015rvir_list.append(stars4_SFR_015rvir)
            stars4_SFR_020rvir_list.append(stars4_SFR_020rvir)
            
            stars0_SSFR_010rvir_list.append(stars0_SSFR_010rvir)
            stars0_SSFR_015rvir_list.append(stars0_SSFR_015rvir)
            stars0_SSFR_020rvir_list.append(stars0_SSFR_020rvir)
            stars1_SSFR_010rvir_list.append(stars1_SSFR_010rvir)
            stars1_SSFR_015rvir_list.append(stars1_SSFR_015rvir)
            stars1_SSFR_020rvir_list.append(stars1_SSFR_020rvir)
            stars2_SSFR_010rvir_list.append(stars2_SSFR_010rvir)
            stars2_SSFR_015rvir_list.append(stars2_SSFR_015rvir)
            stars2_SSFR_020rvir_list.append(stars2_SSFR_020rvir)
            stars3_SSFR_010rvir_list.append(stars3_SSFR_010rvir)
            stars3_SSFR_015rvir_list.append(stars3_SSFR_015rvir)
            stars3_SSFR_020rvir_list.append(stars3_SSFR_020rvir)
            stars4_SSFR_010rvir_list.append(stars4_SSFR_010rvir)
            stars4_SSFR_015rvir_list.append(stars4_SSFR_015rvir)
            stars4_SSFR_020rvir_list.append(stars4_SSFR_020rvir)
        
            rvir1_mass_contamination_list.append(rvir1_mass_contamination)
            rvir_number_contamination_list.append(rvir1_number_contamination)
            rvir1_mass_contamination_fraction_list.append(rvir1_mass_contamination_fraction)
            rvir1_number_contamination_fraction_list.append(rvir1_number_contamination_fraction)
            rvir1_larger_darkmatter_present_list.append(rvir1_larger_darkmatter_present)
            rvir010_number_contamination_list.append(rvir010_number_contamination)
            
            #what we added for 100 Myr
            
            stars100Myr_SFR_010rvir_list.append(stars100Myr_SFR_010rvir)  
            stars100Myr_SFR_015rvir_list.append(stars100Myr_SFR_015rvir)
            stars100Myr_SFR_020rvir_list.append(stars100Myr_SFR_020rvir) 
            
            stars100Myr_SSFR_010rvir_list.append(stars100Myr_SSFR_010rvir)
            stars100Myr_SSFR_015rvir_list.append(stars100Myr_SSFR_015rvir) 
            stars100Myr_SSFR_020rvir_list.append(stars100Myr_SSFR_020rvir)
            
            stars60Myr_SFR_010rvir_list.append(stars60Myr_SFR_010rvir)  
            stars60Myr_SFR_015rvir_list.append(stars60Myr_SFR_015rvir)
            stars60Myr_SFR_020rvir_list.append(stars60Myr_SFR_020rvir) 
            
            stars60Myr_SSFR_010rvir_list.append(stars60Myr_SSFR_010rvir)
            stars60Myr_SSFR_015rvir_list.append(stars60Myr_SSFR_015rvir) 
            stars60Myr_SSFR_020rvir_list.append(stars60Myr_SSFR_020rvir)
            
            #the tacchella numbers 
            
            Tacchella5Myr_SFR_010_final_average_list.append(Tacchella5Myr_SFR_010_final_average)
            Tacchella5Myr_SFR_010_final_median_list.append(Tacchella5Myr_SFR_010_final_median)
            Tacchella5Myr_SSFR_010_final_average_list.append(Tacchella5Myr_SSFR_010_final_average)
            Tacchella5Myr_SSFR_010_final_median_list.append(Tacchella5Myr_SSFR_010_final_median)

            Tacchella5Myr_SFR_015_final_average_list.append(Tacchella5Myr_SFR_015_final_average)
            Tacchella5Myr_SFR_015_final_median_list.append(Tacchella5Myr_SFR_015_final_median)
            Tacchella5Myr_SSFR_015_final_average_list.append(Tacchella5Myr_SSFR_015_final_average)
            Tacchella5Myr_SSFR_015_final_median_list.append(Tacchella5Myr_SSFR_015_final_median)

            Tacchella5Myr_SFR_020_final_average_list.append(Tacchella5Myr_SFR_020_final_average)
            Tacchella5Myr_SFR_020_final_median_list.append(Tacchella5Myr_SFR_020_final_median)
            Tacchella5Myr_SSFR_020_final_average_list.append(Tacchella5Myr_SSFR_020_final_average)
            Tacchella5Myr_SSFR_020_final_median_list.append(Tacchella5Myr_SSFR_020_final_median)
            
            Tacchella60Myr_SFR_010_final_average_list.append(Tacchella60Myr_SFR_010_final_average)
            Tacchella60Myr_SFR_010_final_median_list.append(Tacchella60Myr_SFR_010_final_median)
            Tacchella60Myr_SSFR_010_final_average_list.append(Tacchella60Myr_SSFR_010_final_average)
            Tacchella60Myr_SSFR_010_final_median_list.append(Tacchella60Myr_SSFR_010_final_median)

            Tacchella60Myr_SFR_015_final_average_list.append(Tacchella60Myr_SFR_015_final_average)
            Tacchella60Myr_SFR_015_final_median_list.append(Tacchella60Myr_SFR_015_final_median)
            Tacchella60Myr_SSFR_015_final_average_list.append(Tacchella60Myr_SSFR_015_final_average)
            Tacchella60Myr_SSFR_015_final_median_list.append(Tacchella60Myr_SSFR_015_final_median)

            Tacchella60Myr_SFR_020_final_average_list.append(Tacchella60Myr_SFR_020_final_average)
            Tacchella60Myr_SFR_020_final_median_list.append(Tacchella60Myr_SFR_020_final_median)
            Tacchella60Myr_SSFR_020_final_average_list.append(Tacchella60Myr_SSFR_020_final_average)
            Tacchella60Myr_SSFR_020_final_median_list.append(Tacchella60Myr_SSFR_020_final_median)
            
            Tacchella100Myr_SFR_010_final_average_list.append(Tacchella100Myr_SFR_010_final_average)
            Tacchella100Myr_SFR_010_final_median_list.append(Tacchella100Myr_SFR_010_final_median)
            Tacchella100Myr_SSFR_010_final_average_list.append(Tacchella100Myr_SSFR_010_final_average)
            Tacchella100Myr_SSFR_010_final_median_list.append(Tacchella100Myr_SSFR_010_final_median)

            Tacchella100Myr_SFR_015_final_average_list.append(Tacchella100Myr_SFR_015_final_average)
            Tacchella100Myr_SFR_015_final_median_list.append(Tacchella100Myr_SFR_015_final_median)
            Tacchella100Myr_SSFR_015_final_average_list.append(Tacchella100Myr_SSFR_015_final_average)
            Tacchella100Myr_SSFR_015_final_median_list.append(Tacchella100Myr_SSFR_015_final_median)

            Tacchella100Myr_SFR_020_final_average_list.append(Tacchella100Myr_SFR_020_final_average)
            Tacchella100Myr_SFR_020_final_median_list.append(Tacchella100Myr_SFR_020_final_median)
            Tacchella100Myr_SSFR_020_final_average_list.append(Tacchella100Myr_SSFR_020_final_average)
            Tacchella100Myr_SSFR_020_final_median_list.append(Tacchella100Myr_SSFR_020_final_median)
            
        #now save these to an ascii file
        file_name = '%s/SFRandcontamination%s.ascii' % (out_dir, current_scale)
        
        names_list = ('Halo_Id[0]', 'Rvir(Kpc)[1]', 'Mvir(Msun)[2]',
                     'stars0_SFR_0.1Rvir[3]', 'stars1_SFR_0.1Rvir[4]', 'stars2_SFR_0.1Rvir[5]', 'stars3_SFR_0.1Rvir[6]', 'stars4_SFR_0.1Rvir[7]',
                     'stars0_SSFR_0.1Rvir[8]', 'stars1_SSFR_0.1Rvir[9]', 'stars2_SSFR_0.1Rvir[10]', 'stars3_SSFR_0.1Rvir[11]', 'stars4_SSFR_0.1Rvir[12]',
                     'stars0_SFR_0.15Rvir[13]', 'stars1_SFR_0.15Rvir[14]', 'stars2_SFR_0.15Rvir[15]', 'stars3_SFR_0.15Rvir[16]', 'stars4_SFR_0.15Rvir[17]',
                     'stars0_SSFR_0.1Rvir[18]', 'stars1_SSFR_0.15Rvir[19]', 'stars2_SSFR_0.15Rvir[20]', 'stars3_SSFR_0.15Rvir[21]', 'stars4_SSFR_0.15Rvir[22]',
                     'stars0_SFR_0.2Rvir[23]', 'stars1_SFR_0.2Rvir[24]', 'stars2_SFR_0.2Rvir[25]', 'stars3_SFR_0.2Rvir[26]', 'stars4_SFR_0.2Rvir[27]',
                     'stars0_SSFR_0.2Rvir[28]', 'stars1_SSFR_0.2Rvir[29]', 'stars2_SSFR_0.2Rvir[30]', 'stars3_SSFR_0.2Rvir[31]', 'stars4_SSFR_0.2Rvir[32]',
                     'DM_Contamination_1Rvir(Msun)[33]','DM_Contamination_1Rvir_MassFraction[34]', 'DM_Contamination_1Rvir_Number_Particles[35]', 
                     'DM_Contamination_1Rvir_Number_Particles_Fraction[36]', 'DM_Contamination_LargerParticlesPresent(Bool)[37]', 'DM_Contamination_0.1Rvir_Number_Particles[38]',
                     '60Myr_SFR_0.1Rvir[39]', '60Myr_SFR_0.15Rvir[60]', '60Myr_SFR_0.2Rvir[41]', '60Myr_SSFR_0.1Rvir[42]', '60Myr_SSFR_0.15Rvir[43]', '60Myr_SSFR_0.2Rvir[44]', 
                     '100Myr_SFR_0.1Rvir[45]', '60Myr_SFR_0.15Rvir[46]', '100Myr_SFR_0.2Rvir[47]', '100Myr_SSFR_0.1Rvir[48]', '100Myr_SSFR_0.15Rvir[49]', '100Myr_SSFR_0.2Rvir[50]', 
                     'Tacchella_5Myr_SFR_010_AVG[51]', 'Tacchella_5Myr_SFR_010_MED[52]', 'Tacchella_5Myr_SSFR_010_AVG[53]', 'Tacchella_5Myr_SSFR_010_MED[54]', 
                     'Tacchella_5Myr_SFR_015_AVG[55]', 'Tacchella_5Myr_SFR_015_MED[56]', 'Tacchella_5Myr_SSFR_015_AVG[57]', 'Tacchella_5Myr_SSFR_015_MED[58]', 
                     'Tacchella_5Myr_SFR_020_AVG[59]', 'Tacchella_5Myr_SFR_020_MED[5]', 'Tacchella_5Myr_SSFR_020_AVG[61]', 'Tacchella_5Myr_SSFR_020_MED[62]', 
                     'Tacchella_60Myr_SFR_010_AVG[63]', 'Tacchella_60Myr_SFR_010_MED[64]', 'Tacchella_60Myr_SSFR_010_AVG[65]', 'Tacchella_60Myr_SSFR_010_MED[66]', 
                     'Tacchella_60Myr_SFR_015_AVG[67]', 'Tacchella_60Myr_SFR_015_MED[68]', 'Tacchella_60Myr_SSFR_015_AVG[69]', 'Tacchella_60Myr_SSFR_015_MED[70]', 
                     'Tacchella_60Myr_SFR_020_AVG[71]', 'Tacchella_60Myr_SFR_020_MED[72]', 'Tacchella_60Myr_SSFR_020_AVG[73]', 'Tacchella_60Myr_SSFR_020_MED[74]', 
                     'Tacchella_100Myr_SFR_010_AVG[75]', 'Tacchella_100Myr_SFR_010_MED[76]', 'Tacchella_100Myr_SSFR_010_AVG[77]', 'Tacchella_100Myr_SSFR_010_MED[78]', 
                     'Tacchella_100Myr_SFR_015_AVG[79]', 'Tacchella_100Myr_SFR_015_MED[80]', 'Tacchella_100Myr_SSFR_015_AVG[81]', 'Tacchella_100Myr_SSFR_015_MED[82]', 
                     'Tacchella_100Myr_SFR_020_AVG[83]', 'Tacchella_100Myr_SFR_020_MED[84]', 'Tacchella_100Myr_SSFR_020_AVG[85]', 'Tacchella_100Myr_SSFR_020_MED[86]')
            
        data = Table([halo_id_list, rvir_list, mvir_list, stars0_SFR_010rvir_list, stars1_SFR_010rvir_list, 
                      stars2_SFR_010rvir_list, stars3_SFR_010rvir_list, stars4_SFR_010rvir_list, 
                      stars0_SSFR_010rvir_list, stars1_SSFR_010rvir_list, stars2_SSFR_010rvir_list, 
                      stars3_SSFR_010rvir_list, stars4_SSFR_010rvir_list, stars0_SFR_015rvir_list, 
                      stars1_SFR_015rvir_list, stars2_SFR_015rvir_list, stars3_SFR_015rvir_list, 
                      stars4_SFR_015rvir_list, stars0_SSFR_015rvir_list, stars1_SSFR_015rvir_list, 
                      stars2_SSFR_015rvir_list, stars3_SSFR_015rvir_list, stars4_SSFR_015rvir_list, 
                      stars0_SFR_020rvir_list, stars1_SFR_020rvir_list, stars2_SFR_020rvir_list, 
                      stars3_SFR_020rvir_list, stars4_SFR_020rvir_list, stars0_SSFR_020rvir_list, 
                      stars1_SSFR_020rvir_list, stars2_SSFR_020rvir_list, stars3_SSFR_020rvir_list, 
                      stars4_SSFR_020rvir_list, rvir1_mass_contamination_list, rvir1_mass_contamination_fraction_list, 
                      rvir_number_contamination_list, rvir1_number_contamination_fraction_list, 
                      rvir1_larger_darkmatter_present_list, rvir010_number_contamination_list, 
                      stars60Myr_SFR_010rvir_list, stars60Myr_SFR_015rvir_list, stars60Myr_SFR_020rvir_list, 
                      stars60Myr_SSFR_010rvir_list, stars60Myr_SSFR_015rvir_list, stars60Myr_SSFR_020rvir_list, 
                      stars100Myr_SFR_010rvir_list, stars100Myr_SFR_015rvir_list, stars100Myr_SFR_020rvir_list, 
                      stars100Myr_SSFR_010rvir_list, stars100Myr_SSFR_015rvir_list, stars100Myr_SSFR_020rvir_list, 
                      Tacchella5Myr_SFR_010_final_average_list, Tacchella5Myr_SFR_010_final_median_list, 
                      Tacchella5Myr_SSFR_010_final_average_list, Tacchella5Myr_SSFR_010_final_median_list, 
                      Tacchella5Myr_SFR_015_final_average_list, Tacchella5Myr_SFR_015_final_median_list,
                      Tacchella5Myr_SSFR_015_final_average_list, Tacchella5Myr_SSFR_015_final_median_list, 
                      Tacchella5Myr_SFR_020_final_average_list, Tacchella5Myr_SFR_020_final_median_list, 
                      Tacchella5Myr_SSFR_020_final_average_list, Tacchella5Myr_SSFR_020_final_median_list, 
                      Tacchella60Myr_SFR_010_final_average_list, Tacchella60Myr_SFR_010_final_median_list, 
                      Tacchella60Myr_SSFR_010_final_average_list, Tacchella60Myr_SSFR_010_final_median_list, 
                      Tacchella60Myr_SFR_015_final_average_list, Tacchella60Myr_SFR_015_final_median_list, 
                      Tacchella60Myr_SSFR_015_final_average_list, Tacchella60Myr_SSFR_015_final_median_list, 
                      Tacchella60Myr_SFR_020_final_average_list, Tacchella60Myr_SFR_020_final_median_list, 
                      Tacchella60Myr_SSFR_020_final_average_list, Tacchella60Myr_SSFR_020_final_median_list, 
                      Tacchella100Myr_SFR_010_final_average_list, Tacchella100Myr_SFR_010_final_median_list, 
                      Tacchella100Myr_SSFR_010_final_average_list, Tacchella100Myr_SSFR_010_final_median_list, 
                      Tacchella100Myr_SFR_015_final_average_list, Tacchella100Myr_SFR_015_final_median_list, 
                      Tacchella100Myr_SSFR_015_final_average_list, Tacchella100Myr_SSFR_015_final_median_list, 
                      Tacchella100Myr_SFR_020_final_average_list, Tacchella100Myr_SFR_020_final_median_list, 
                      Tacchella100Myr_SSFR_020_final_average_list, Tacchella100Myr_SSFR_020_final_median_list], names = names_list)
        ascii.write(data, output=file_name, overwrite=True)   

