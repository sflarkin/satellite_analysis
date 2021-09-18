def nircatalogreader(GEN, VELA_number):
    if GEN == 3:
        input_dir = '/Users/user1/documents/GEN3/nircatalogs/galaxy_catalogs'
        catalog_dir = '{}/VELA_v2_{}/galaxy_catalogue/'.format(input_dir, VELA_number)
    if GEN == 6:
        input_dir = '/Users/user1/documents/GEN6/nircatalogs/gen6_catalogs'
        catalog_dir = '{}/V{}/galaxy_catalogue/'.format(input_dir, VELA_number)
    
    
    #I am only interested in a few of the files of Nir's catalogs, so I put the ones I want in the following list
    #if you want more of them, add them here, and they will be added in order
    
    wanted_files = ['Nir_simplified_disc_cat.txt', 'Nir_halo_cat.txt', 'Nir_015_Rvir_cat.txt', 'Mstar.txt', 'Nir_disc_cat.txt']
    
    #this is a variable that checks the number of snaps, to check if each file read has the same number of snapshots as a check to 
    #make sure nothing went wrong
    number_snaps = False
    
    #now I read the files I am interested in, and store the info I am interested in
    
    for file in wanted_files:
    
        readfile = open('{}{}'.format(catalog_dir, file))
        lines = readfile.readlines()
    
        catalog = []
        for line in lines:
            catalog.append(line.split())
        
        #if the list of lists to store the info has yet to be formed, create it and set the number_snaps. If it has been created
        #make sure the number_snaps is correct for the files
        
        if number_snaps == False:
            number_snaps = int(catalog[0][0])
            nir_data = [[] for x in range(number_snaps)]
            indices = [x for x in range(number_snaps)]
        else:
            if number_snaps != int(catalog[0][0]):
                print('Number of snaps does not match for file {}'.format(file))
                                   
        #remove the first line of the catalog, that tells you how many snaps there are
        del catalog[0]
    
        for index in indices:
            nir_data[index].append(catalog[index])
        readfile.close()
    return(nir_data)
    
    
    
        

                             
                    