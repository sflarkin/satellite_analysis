import glob
import numpy as np
from operator import itemgetter

def consistent_catalog_reader(input_dir, add_all=False, remove_subhalos='False', halo_mass=1e+08, halo_number=1):
    global halo_data_sorted, halo_data_all, snapshot_index, consistent_file_index, consistent_index
    
    all_files = glob.glob(input_dir + '/*.list')
    consistent_file_index_list = []
    for file in all_files:
        period = [pos for pos, char in enumerate(file) if char == '.']
        digits = file[period[-2]+1:period[-2]+4]
        consistent_file_index_list.append(digits)  

    consistent_file_index = list(set(consistent_file_index_list))
    consistent_file_index.sort()
    
    
    snapshot_index = [x for x in range(len(consistent_file_index))]
    
    print('Consistent File Indices:', consistent_file_index)
    print('')
    
    print('Using Indices:', snapshot_index)
    print('')
    
    #this is the index of what each column is, will be created with the first list
    consistent_index = []

    #create the list of lists that the halo data will be added to

    halo_data_all = [[] for x in range(len(snapshot_index))]

    #create the list of lists for the id's of the halo(s) with the most particles and highest mvir
    halo_data_sorted = [[] for x in range(len(snapshot_index))]
    
    for index in snapshot_index:
    
        #print('')
        #print('Collecting Halo Data for snapshot index', index)
        glob_files = glob.glob(input_dir + '/*hlist_0.' + consistent_file_index[index] + '*.list')

        #call the empty list from the halo_data_lists that we want to add to
    
        halo_list = []
    
        #this opens up each individual snapshot file, deletes the info lines, and adds all of the individual halo
        #data to a listname which corresponds to halos + snapshot_index(halo5 for snapshot 5, for example),
        #iterating over all of the snapshot files for each snapshot
    
        for file in glob_files:
            readfile = open(file)
            lines = readfile.readlines()
            catalog = []
            above_halo_mass = []
            #set the consistent_index if it hasn't already been set
            if consistent_index == []:
                consistent_index = lines[0]
            for line in lines:
                catalog.append(line.split())
            del catalog[0:65]
            if catalog != []:
                #part that derermines if mass is a discriminator or no discriminator
                if add_all == True:
                    halo_list = halo_list + catalog
                if add_all == False:
                    for lists in catalog:
                        if float(lists[10]) >= halo_mass:
                            if remove_subhalos == 'False':
                                above_halo_mass.append(lists)
                            if remove_subhalos == 'True':
                                if float(lists[5]) == -1:
                                    above_halo_mass.append(lists)
                    halo_list = halo_list + above_halo_mass
            readfile.close()
        print('Number of Halos found for snapshot', index, ':', len(halo_list))
        halo_data_all[index] = halo_data_all[index] + halo_list
        
        #now we sort the halo data by Mvir to be stored in 
   

        #halo_data_all[index]   
        mvir_sort = sorted(halo_data_all[index], key=lambda x: float(x[10]), reverse=True)
         
        #now we get the halo data from these largest halos ids and add the info to halo_data_num_p_mvir
        halo_data_sorted[index] = mvir_sort
