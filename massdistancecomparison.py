import numpy as np
import argparse
import matplotlib.pyplot as plt
from satellite_analysis import rockstarcatalogreader


def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('out_file')
    args = vars(parser.parse_args())
    return args

args = parse()
input_dir = args['input_dir']

rockstarcatalogreader.rockstar_catalog_reader(input_dir)

wanted_a = ['170', '200', '250', '330', '400', '500'] #needs to be *1000 since the rockstar_files are like that
scale_indexes = [pos for pos, char in enumerate(rockstarcatalogreader.rockstar_file_index) if char in wanted_a]
plt.figure(figsize=(20,30))
count = 1
for position in scale_indexes:
    
    largest_halo = rockstarcatalgoreader.halo_data_num_p_mvir[position][1][0]
    z = ((1/(float(rockstarcatalogreader.rockstar_file_index[position])/1000)) - 1 )
    z_str = str(z)
    distances = []
    masses = []
    for halo in rockstarcatalogreader.halo_data_all[position]:
        x = (float(largest_halo[8]) - float(halo[8]))
        y = (float(largest_halo[9]) - float(halo[9]))
        z = (float(largest_halo[10]) - float(halo[10]))
        mass = float(halo[2])
        distances.append(((x**2 + y**2 + z**2)**(1/2)))
        masses.append(mass)
    masses_sort = sorted(masses)
    distances_sort = sorted(distances)
    title = 'Z = ' + z_str
    
    plt.subplot(6,3,count)
    plt.xlabel('Mass in MSun/h')
    plt.ylabel('Objects in Mass Bin')
    plt.xlim(10**(7.8), 10**(12.5))
    plt.xscale('log', nonposy='clip')
    plt.yscale('log', nonposy='clip')
    bins = np.logspace(np.log10(10e7),np.log10(10e12), 50)
    n, bins,patches = plt.hist(masses, bins=bins, facecolor='g', alpha = 0.75)
    count = count + 1
    
    num_objects = len(masses)
    text = 'Total Objects:' + str(num_objects)
    plt.subplot(6,3,count)
    plt.xlim(-.025,1)
    plt.ylim(10**(7.8), 10**(12.5))
    plt.yscale('log')
    plt.ylabel('Mass in MSun/h')
    plt.xlabel('Distance in Mpc/h')
    plt.scatter(distances, masses, marker='.', s =.5)
    plt.title(title)
    plt.annotate(text, (.6, 10**12))
    count = count + 1
    
    plt.subplot(6,3,count)
    plt.xlabel('Distance in Mpc/h')
    plt.xlim(-.025,1)
    plt.ylim(0,140)
    plt.ylabel('Objects in Distance Bin')
    n, bins, patches = plt.hist(distances, 50, facecolor='tab:orange', alpha=0.75)
    count = count + 1

plt.savefig('%s/mass&distancecomparison.jpg' % input_dir)  
    
