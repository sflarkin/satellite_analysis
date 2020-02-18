import yt
import numpy as np
import glob
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

file_names = glob.glob('/la1/dceverin/VELA_v2.0/VELA*/timing.log')

list_of_wanted_snaps = [i for i in np.arange(.100, .5025, .0025)]

for file in file_names:
    readfile = open(file)
    lines = readfile.readlines()
    catalog = []
    for line in lines:
        catalog.append(line.split())
    del catalog[0:2]
    print(catalog[0])
    
    readfile.close()
    
    snaps_for_analysis = []
    ids = []
    scales = []
    for wanted_snap in list_of_wanted_snaps:
        time_difference = 10000
        clocest_snap = 0
        for sim_snaps in catalog:
            difference = abs(float(sim_snaps[2]) - wanted_snap)
            if difference < time_difference:
                time_difference = difference
                clocest_snap = (sim_snaps[0], sim_snaps[2])
        ids.append(clocest_snap[0])
        scales.append(clocest_snap[1])
    
    #now make the info to write to the file
    id_string_formatted = ''
    for id_ in ids:
        id_formatted = ' *{}*'.format(('0' * (5 -len(str(id_)))) + str(id_))
        id_string_formatted = id_string_formatted + id_formatted
    shiftc_command = 'shiftc -r {} pfe:/nobackupp2/sflarkin/GEN3/VELA{}'.format(id_string_formatted, file[-13:-11])
    
    comment = {'comments':[shiftc_command]} 
    
    file_name = '/u/sflarkin/VELA{}_snaps.ascii'.format(file[-13:-11])
    
    #create the ascii table
    data = Table([ids, scale], names=['ID', 'scale'], meta=comment)

    ascii.write(data, output=file_name, comment='#')                    
        