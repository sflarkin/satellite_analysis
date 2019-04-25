import glob
import os, sys, argparse
import re
import subprocess
import time

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('--output_dir', nargs='?', default='subhalos')
    parser.add_argument('--rockstar_dir', nargs='?', default='/Users/user1/rockstar')
    args = vars(parser.parse_args())
    return args

args = parse()
#input_dir need to be of form /path/to/baserockstar_ascii
input_dir = args['input_dir']
print('input_dir: %s' % input_dir)

output_dir = args['output_dir']
output_directory = '%s/%s' % (input_dir, output_dir)
print('output_directory: %s' % output_directory)

rockstar_dir = args['rockstar_dir']
print('rockstar_dir: %s' % rockstar_dir)

#Check to see if the output directory exists. If it does not, create it.
if not os.path.exists(output_directory):
    print('Creating output directory')
    os.makedirs(output_directory)

#Find the filenames for all of the rockstar ascii files.
rockstar_files = glob.glob('%s/*.ascii' % input_dir)

for file in rockstar_files:
    #Generate the box size from the a of the rockstar file. 
    search_a = re.search('halos_(...).', file)
    box_size = int(search_a.group(1))/1000*40
    
    #Make the string for the output location for the rockstar/util/find_parents to send to
    slash = [pos for pos, char in enumerate(file) if char == '/']
    halos_ = file[slash[-1]+1:]
    output_file = '>%s/sub%s' % (output_directory, halos_)
    print('output_file:', output_file )   
    #Create the command that will be passed to subprocess.run
    command = '%s/util/find_parents %s %s %s' % (rockstar_dir, file, str(box_size), output_file)
    #used shell=True since using command = ['%s/util/find_parents' % rockstar_dir, file, str(box_size), output_file)]
    #with shell=False did not save the output to the output_file, it just returned it to the command line
    #
    #look into why that is
    subprocess.run(command, shell=True)
    time.sleep(4)
    
print('Successfully completed subhalo catalog for %s' % input_dir)
