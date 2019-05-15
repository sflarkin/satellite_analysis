import glob
import numpy as np
from operator import itemgetter

def consistent_catalog_reader(input_dir, add_all=False, subhalos='False', halo_mass=1e+08, halo_number=20):
    global halo_data_largest, halo_data_all, snapshot_index, consistent_file_index, consistent_index
    
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
    halo_data_largest = [[] for x in range(len(snapshot_index))]
    
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
                            if subhalos == 'False':
                                above_halo_mass.append(lists)
                            if subhalos == 'True':
                                if float(lists[5]) == -1:
                                    above_halo_mass.append(lists)
                    halo_list = halo_list + above_halo_mass
            readfile.close()
        print('Number of Halos found for snapshot', index, ':', len(halo_list))
        halo_data_all[index] = halo_data_all[index] + halo_list
        
        #now we will find the halo(s) with the highest num_p and mvir from the halo_list
    
        #empty the lists
        id_mvir_rvir = []
        mvir_sort, rvir_sort = [], []
        mvir_ids, rvir_ids = [], []
        mvir_list, rvir_list = [], []
        mvir_list_sort, rvir_list_sort = [], []

        #halo_data_all[index]
        for lines in halo_data_all[index]:
            id_mvir_rvir.append([float(lines[1]), float(lines[10]), float(lines[11])])     
        mvir_sort = sorted(id_mvir_rvir, key=itemgetter(1), reverse=True)
        rvir_sort = sorted(id_mvir_rvir, key=itemgetter(2), reverse=True)
                        
        if mvir_sort != []:
            if len(mvir_sort) <= halo_number:
                mvir_ids = [mvir_sort[y][0] for y in range(len(mvir_sort))]
            else:
                mvir_ids = [mvir_sort[y][0] for y in range(halo_number)]
            
            for lines in halo_data_all[index]:
                if float(lines[1]) in mvir_ids:
                    mvir_list.append(lines)
            for ids in mvir_ids:
                for lines in mvir_list:
                    if str(int(ids)) == lines[1]:
                        mvir_list_sort.append(lines)
        
        if rvir_sort != []:
            if len(rvir_sort) <= halo_number:
                rvir_ids = [rvir_sort[y][0] for y in range(len(rvir_sort))]
            else:
                rvir_ids = [rvir_sort[y][0] for y in range(halo_number)]
            
            for lines in halo_data_all[index]:
                if float(lines[1]) in rvir_ids:
                    rvir_list.append(lines)
            for ids in rvir_ids:
                for lines in rvir_list:
                    if str(int(ids)) == lines[1]:
                        rvir_list_sort.append(lines)
                        
        #now we get the halo data from these largest halos ids and add the info to halo_data_num_p_mvir
        halo_data_largest[index] = [mvir_list_sort, rvir_list_sort]
