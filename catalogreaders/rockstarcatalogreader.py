import glob
import numpy as np
from operator import itemgetter

def rockstar_catalog_reader(input_dir, add_all=False, subhalos='False', halo_mass=1e+08, halo_number=20):
    global halo_data_largest, halo_data_all, snapshot_index, rockstar_file_index, rockstar_index
    
    all_files = glob.glob(input_dir + '/*.ascii')
    rockstar_file_index_list = []
    for file in all_files:
        underscore = [pos for pos, char in enumerate(file) if char == '_']
        period = [pos for pos, char in enumerate(file) if char == '.']
        digits = file[underscore[-1]+1:period[-2]]
        rockstar_file_index_list.append(digits)  

    rockstar_file_index = list(set(rockstar_file_index_list))
    rockstar_file_index.sort()
    
    snapshot_index = [x for x in range(len(rockstar_file_index))]
    
    print('Rockstar File Indices:', rockstar_file_index)
    print('')
    
    print('Using Indices:', snapshot_index)
    print('')
    
    #this is the index of what each column is, will be created with the first list
    rockstar_index = []

    #create the list of lists that the halo data will be added to

    halo_data_all = [[] for x in range(len(snapshot_index))]

    #create the list of lists for the id's of the halo(s) with the most particles and highest mvir
    halo_data_largest = [[] for x in range(len(snapshot_index))]
    
    for index in snapshot_index:
    
        #print('')
        #print('Collecting Halo Data for snapshot index', index)
        glob_files = glob.glob(input_dir + '/*halos_' + rockstar_file_index[index] + '.*.ascii')

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
            #set the rockstar_index if it hasn't already been set
            if rockstar_index == []:
                rockstar_index = lines[0]
            for line in lines:
                catalog.append(line.split())
            del catalog[0:20]
            if catalog != []:
                #part that derermines if mass is a discriminator or no discriminator
                if add_all == True:
                    halo_list = halo_list + catalog
                if add_all == False:
                    for lists in catalog:
                        if float(lists[2]) >= halo_mass:
                            if subhalos == 'False':
                                above_halo_mass.append(lists)
                            if subhalos == 'True':
                                if float(lists[-1]) == -1:
                                    above_halo_mass.append(lists)
                    halo_list = halo_list + above_halo_mass
            readfile.close()
        print('Number of Halos found for snapshot', index, ':', len(halo_list))
        halo_data_all[index] = halo_data_all[index] + halo_list
        
        #now we will find the halo(s) with the highest num_p and mvir from the halo_list
    
        #empty the lists
        id_num_p_mvir = []
        num_p_sort, mvir_sort, rvir_sort = [], [], []
        num_p_ids, mvir_ids, rvir_ids = [], [], []
        num_p_list, mvir_list, rvir_list = [], [], []
        num_p_list_sort, mvir_list_sort, rvir_list_sort = [], [], []

        #halo_data_all[index]
        for lines in halo_data_all[index]:
            id_num_p_mvir.append([float(lines[0]), float(lines[1]), float(lines[2]), float(lines[4])])     
        num_p_sort = sorted(id_num_p_mvir, key=itemgetter(1), reverse=True)
        mvir_sort = sorted(id_num_p_mvir, key=itemgetter(2), reverse=True)
        rvir_sort = sorted(id_num_p_mvir, key=itemgetter(3), reverse=True)
        
        if num_p_sort != []:
            if len(num_p_sort) <= halo_number:
                num_p_ids = [num_p_sort[y][0] for y in range(len(num_p_sort))]
            else:
                num_p_ids = [num_p_sort[y][0] for y in range(halo_number)]
            
            for lines in halo_data_all[index]:
                if float(lines[0]) in num_p_ids:
                    num_p_list.append(lines)
            for ids in num_p_ids:
                for lines in num_p_list:
                    if str(int(ids)) == lines[0]:
                        num_p_list_sort.append(lines)
                        
        if mvir_sort != []:
            if len(mvir_sort) <= halo_number:
                mvir_ids = [mvir_sort[y][0] for y in range(len(mvir_sort))]
            else:
                mvir_ids = [mvir_sort[y][0] for y in range(halo_number)]
            
            for lines in halo_data_all[index]:
                if float(lines[0]) in mvir_ids:
                    mvir_list.append(lines)
            for ids in mvir_ids:
                for lines in mvir_list:
                    if str(int(ids)) == lines[0]:
                        mvir_list_sort.append(lines)
        
        if rvir_sort != []:
            if len(rvir_sort) <= halo_number:
                rvir_ids = [rvir_sort[y][0] for y in range(len(rvir_sort))]
            else:
                rvir_ids = [rvir_sort[y][0] for y in range(halo_number)]
            
            for lines in halo_data_all[index]:
                if float(lines[0]) in rvir_ids:
                    rvir_list.append(lines)
            for ids in rvir_ids:
                for lines in rvir_list:
                    if str(int(ids)) == lines[0]:
                        rvir_list_sort.append(lines)
                        
        #now we get the halo data from these largest halos ids and add the info to halo_data_num_p_mvir
        halo_data_largest[index] = [num_p_list_sort, mvir_list_sort, rvir_list_sort]
        