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

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('proper_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--remove_subhalos', nargs='?', default='True')
    parser.add_argument('--halo_mass', nargs='?', default=10e+09)
    args = vars(parser.parse_args())
    return args


args = parse()

consistent_dir = args['input_dir']
proper_dir = args['proper_dir']
VELA_dir = args['VELA_dir']
out_dir = args['out_dir']
remove_subhalos = args['remove_subhalos']
halo_mass = args['halo_mass']

completed = glob.glob('%s/o2abundance*.ascii' % (out_dir))


#Check if the output directory exists and if it does not, create it.
if not os.path.exists(out_dir):
    print('Creating output directory')
    os.makedirs(out_dir)

VELA_snaps = glob.glob(VELA_dir + '/10MpcBox*')
VELA_snaps.sort()
print(VELA_snaps)

consistent.consistent_catalog_reader(consistent_dir, remove_subhalos='True', halo_mass=10e9)

def _massSNII(field, data):
    return data['MetalDensitySNII'] * data['cell_volume']
yt.add_field(('art', 'MetalMassSNII'), function=_massSNII, units='g')

def _massSNIa(field, data):
    return data['MetalDensitySNIa'] * data['cell_volume']
yt.add_field(('art', 'MetalMassSNIa'), function=_massSNIa, units='g')

for snap in reversed(VELA_snaps):
    ds = yt.load(snap)
    current_scale = str(ds.scale_factor)[2:5]
    
    current_title = '%s/o2abundance%s.ascii' % (out_dir, current_scale)
    if current_title in completed:
        continue 
        
    print('Trying to find matching Consistent File for scale: {}'.format(current_scale))
    
    position = [pos for pos, loc in enumerate(consistent.consistent_file_index) if loc == current_scale]

    if position == [] or len(position) > 1:
        print('Could not find corresponding consistent trees data for file {}'.format(snap))
    else:
        
        print('Finding MVir Masses for snap:{}'.format(snap))
        print('Matching Scales Test:{} {}'.format(current_scale, consistent.consistent_file_index[position[0]]))
        
        #make the lists to hold the data to be written to files later
        
        rockstar_id_list = []
        
        proper_rvir_list = []
        
        type2_2percent_list, gas_mass_2percent_list, o2abundance_2percent_list = [], [], [] 
        type2_4percent_list, gas_mass_4percent_list, o2abundance_4percent_list = [], [], []
        type2_6percent_list, gas_mass_6percent_list, o2abundance_6percent_list = [], [], []
        type2_8percent_list, gas_mass_8percent_list, o2abundance_8percent_list = [], [], []
        type2_10percent_list, gas_mass_10percent_list, o2abundance_10percent_list = [], [], []
        type2_15percent_list, gas_mass_15percent_list, o2abundance_15percent_list = [], [], []
        type2_20percent_list, gas_mass_20percent_list, o2abundance_20percent_list = [], [], []
        
        type1a_2percent_list, type1a_4percent_list, type1a_6percent_list = [], [], []
        type1a_8percent_list, type1a_10percent_list, type1a_15percent_list, type1a_20percent_list = [], [], [], []
        
        
        count = 0
        
        halo_data = consistent.halo_data_sorted[position[0]]
        print(halo_data)

        for halo in halo_data:
            
            #set the values to be calculated to None, so that if any values arent found we will know
            
            
            x = ds.quan(float(halo[17]), 'Mpccm/h')
            y = ds.quan(float(halo[18]), 'Mpccm/h')
            z = ds.quan(float(halo[19]), 'Mpccm/h')
            rockstar_id = float(halo[1])
            center = [x, y, z]
            print(center)
            rockstar_rvir = ds.quan(float(halo[11]), 'kpccm/h')
            
            #now we get the rvir from the propermvircalc files
            
            proper_file = '{}/improvedrvir{}.ascii'.format(proper_dir, current_scale)
            
            readfile = open(proper_file)
            lines = readfile.readlines()
            halomass_lines = []
            for line in lines:
                halomass_lines.append(line.split())
            del halomass_lines[0]
            
            for line in halomass_lines:
                if float(line[0]) == rockstar_id:
                    proper_rvir = ds.quan(float(line[3]), 'kpc')
            
            #this is what makes the image
            
            #need to find the rvir for making this
            
            
            
            sp_proj = ds.sphere(center, (float(proper_rvir.in_units('kpc')), 'kpc'))
            
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
                                  ('stars', 'particle_mass'), center=center, width=(float(proper_rvir.in_units('kpc'))*2, "kpc"), data_source=sp_proj)
            b = yt.ParticlePlot(ds, ('stars', 'particle_position_y'), ('stars', 'particle_position_z'),\
                                  ('stars', 'particle_mass'), center=center, width=(float(proper_rvir.in_units('kpc'))*2, "kpc"), data_source=sp_proj)
            c = yt.ParticlePlot(ds, ('stars', 'particle_position_z'), ('stars', 'particle_position_x'),\
                                  ('stars', 'particle_mass'), center=center, width=(float(proper_rvir.in_units('kpc'))*2, "kpc"), data_source=sp_proj)

            d = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_x'), ('darkmatter', 'particle_position_y'),\
                                  ('darkmatter', 'particle_mass'), center=center, width=(float(proper_rvir.in_units('kpc'))*2, "kpc"), data_source=sp_proj)
            e = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_y'), ('darkmatter', 'particle_position_z'),\
                                  ('darkmatter', 'particle_mass'), center=center, width=(float(proper_rvir.in_units('kpc'))*2, "kpc"), data_source=sp_proj)
            f = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_z'), ('darkmatter', 'particle_position_x'),\
                                  ('darkmatter', 'particle_mass'), center=center, width=(float(proper_rvir.in_units('kpc'))*2, "kpc"), data_source=sp_proj)

            a.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            b.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            c.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            d.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            e.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            f.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            
            a.annotate_sphere(center, radius=(proper_rvir.in_units('kpc') * 0.15, 'kpc'), circle_args={'color':'black'})
            b.annotate_sphere(center, radius=(proper_rvir.in_units('kpc') * 0.15, 'kpc'), circle_args={'color':'black'})
            c.annotate_sphere(center, radius=(proper_rvir.in_units('kpc') * 0.15, 'kpc'), circle_args={'color':'black'})
            d.annotate_sphere(center, radius=(rockstar_rvir.in_units('kpc'), 'kpc'), circle_args={'color':'yellow'})
            e.annotate_sphere(center, radius=(rockstar_rvir.in_units('kpc'), 'kpc'), circle_args={'color':'yellow'})
            f.annotate_sphere(center, radius=(rockstar_rvir.in_units('kpc'), 'kpc'), circle_args={'color':'yellow'})
            

            index = 0
            for letter in [a,b,c]:
                plot = letter.plots[('stars', 'particle_mass')]
                plot.figure = fig
                plot.axes = grid[index].axes
                plot.cax = grid.cbar_axes[0]
                letter._setup_plots()
                index = index + 1

            for letter in [d,e,f]:
                plot = letter.plots[('darkmatter', 'particle_mass')]
                plot.figure = fig
                plot.axes = grid[index].axes
                plot.cax = grid.cbar_axes[0]
                letter._setup_plots()
                index = index + 1

            plt.savefig('{}/a{}_{}_id{}.png'.format(out_dir, current_scale, count, rockstar_id), bbox_inches='tight')
            plt.close()

            
            
            #rvirprj.rvirprojections(ds, center, rvir, sp_proj, out_dir, current_scale, rockstar_id, count)
            count = count+1
            #continue
            #now we find the gas fields we are interested in
            #want 2%, 4%, 6%, 8% 10%, 15% and 20% of the virial radius
            
            
            sp_2percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.02, 'kpc'))
            sp_4percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.04, 'kpc'))
            sp_6percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.06, 'kpc'))
            sp_8percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.08, 'kpc'))
            sp_10percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.10, 'kpc'))
            
            sp_15percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.15, 'kpc'))
            sp_20percent = ds.sphere(center, (float(proper_rvir.in_units('kpc'))*.2, 'kpc'))
            
            
            #create the fixed density finder field
            
            
            #now we find the gas mass and type2 supernovae mass for each
            
            gas_mass_2percent, type2_2percent, type1a_2percent = sp_2percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            gas_mass_4percent, type2_4percent, type1a_4percent = sp_4percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            gas_mass_6percent, type2_6percent, type1a_6percent = sp_6percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            gas_mass_8percent, type2_8percent, type1a_8percent = sp_8percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            gas_mass_10percent, type2_10percent, type1a_10percent = sp_10percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            gas_mass_15percent, type2_15percent, type1a_15percent = sp_15percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            gas_mass_20percent, type2_20percent, type1a_20percent = sp_20percent.quantities.total_quantity([('gas', 'cell_mass'), ('art', 'MetalMassSNII'), ('art', 'MetalMassSNIa')])
            
            #the type2 field is a density, so we will have to convert it to mass
            #the *.5 is to only get half which is the part of the metals that are oxygen
            
            mass_o2_2percent = type2_2percent.in_units('Msun') * .5
            mass_o2_4percent = type2_4percent.in_units('Msun') * .5
            mass_o2_6percent = type2_6percent.in_units('Msun') * .5
            mass_o2_8percent = type2_8percent.in_units('Msun') * .5
            mass_o2_10percent = type2_10percent.in_units('Msun') * .5
            mass_o2_15percent = type2_15percent.in_units('Msun') * .5
            mass_o2_20percent = type2_20percent.in_units('Msun') * .5
            
            
            #the * .75 for the gas is because 75% of gas in universe is Hydrogen 
            o2abundance_2percent = 12 + np.log10(mass_o2_2percent / (gas_mass_2percent.in_units('Msun') *.75))
            o2abundance_4percent = 12 + np.log10(mass_o2_4percent / (gas_mass_4percent.in_units('Msun') *.75))
            o2abundance_6percent = 12 + np.log10(mass_o2_6percent / (gas_mass_6percent.in_units('Msun') *.75))
            o2abundance_8percent = 12 + np.log10(mass_o2_8percent / (gas_mass_8percent.in_units('Msun') *.75))
            o2abundance_10percent = 12 + np.log10(mass_o2_10percent / (gas_mass_10percent.in_units('Msun') *.75))
            o2abundance_15percent = 12 + np.log10(mass_o2_15percent / (gas_mass_15percent.in_units('Msun') *.75))
            o2abundance_20percent = 12 + np.log10(mass_o2_20percent / (gas_mass_20percent.in_units('Msun') *.75))
            
            

            #now add the values calculated to the lists for saving to a file
            rockstar_id_list.append(rockstar_id)
        
            proper_rvir_list.append(float(proper_rvir.in_units('kpc')))
            
            type2_2percent_list.append(float(type2_2percent.in_units('Msun')))
            gas_mass_2percent_list.append(float(gas_mass_2percent.in_units('Msun')))
            o2abundance_2percent_list.append(float(o2abundance_2percent))
            
            type2_4percent_list.append(float(type2_4percent.in_units('Msun')))
            gas_mass_4percent_list.append(float(gas_mass_4percent.in_units('Msun')))
            o2abundance_4percent_list.append(float(o2abundance_4percent))
            
            type2_6percent_list.append(float(type2_6percent.in_units('Msun')))
            gas_mass_6percent_list.append(float(gas_mass_6percent.in_units('Msun')))
            o2abundance_6percent_list.append(float(o2abundance_6percent))
            
            type2_8percent_list.append(float(type2_8percent.in_units('Msun')))
            gas_mass_8percent_list.append(float(gas_mass_8percent.in_units('Msun')))
            o2abundance_8percent_list.append(float(o2abundance_8percent))
            
            type2_10percent_list.append(float(type2_10percent.in_units('Msun')))
            gas_mass_10percent_list.append(float(gas_mass_10percent.in_units('Msun')))
            o2abundance_10percent_list.append(float(o2abundance_10percent))
            
            type2_15percent_list.append(float(type2_15percent.in_units('Msun')))
            gas_mass_15percent_list.append(float(gas_mass_15percent.in_units('Msun')))
            o2abundance_15percent_list.append(float(o2abundance_15percent))
            
            type2_20percent_list.append(float(type2_20percent.in_units('Msun')))
            gas_mass_20percent_list.append(float(gas_mass_20percent.in_units('Msun')))
            o2abundance_20percent_list.append(float(o2abundance_20percent))
            
            type1a_2percent_list.append(float(type1a_2percent.in_units('Msun')))
            type1a_4percent_list.append(float(type1a_4percent.in_units('Msun')))
            type1a_6percent_list.append(float(type1a_6percent.in_units('Msun')))
            type1a_8percent_list.append(float(type1a_8percent.in_units('Msun')))
            type1a_10percent_list.append(float(type1a_10percent.in_units('Msun')))
            type1a_15percent_list.append(float(type1a_15percent.in_units('Msun')))
            type1a_20percent_list.append(float(type1a_20percent.in_units('Msun')))
        
            
        #now write the values to a file
        
        ascii_data = [rockstar_id_list, proper_rvir_list, type2_2percent_list, gas_mass_2percent_list, o2abundance_2percent_list, 
                       type2_4percent_list, gas_mass_4percent_list, o2abundance_4percent_list, type2_6percent_list, gas_mass_6percent_list,
                       o2abundance_6percent_list, type2_8percent_list, gas_mass_8percent_list, o2abundance_8percent_list, type2_10percent_list, 
                       gas_mass_10percent_list, o2abundance_10percent_list, type2_15percent_list, gas_mass_15percent_list, o2abundance_15percent_list, 
                       type2_20percent_list, gas_mass_20percent_list, o2abundance_20percent_list, type1a_2percent_list, type1a_4percent_list,
                       type1a_6percent_list, type1a_8percent_list, type1a_10percent_list, type1a_15percent_list, type1a_20percent_list]
        
        ascii_names = ['Rockstar_Id[0]', 'Proper_Rvir[1]', 
                       'SN2_2percent_density[2]', 'Gas_2percent_mass[3]', 'o2_abundance_2percent[4]',
                       'SN2_4percent_density[5]', 'Gas_4percent_mass[6]', 'o2_abundance_4percent[7]', 
                       'SN2_6percent_density[8]', 'Gas_6percent_mass[9]', 'o2_abundance_6percent[10]', 
                       'SN2_8percent_density[11]', 'Gas_8percent_mass[12]', 'o2_abundance_8percent[13]',
                       'SN2_10percent_density[14]', 'Gas_10percent_mass[15]', 'o2_abundance_10percent[16]',
                       'SN2_15percent_density[17]', 'Gas_15percent_mass[18]', 'o2_abundance_15percent[19]',
                       'SN2_20percent_density[20]', 'Gas_20percent_mass[21]', 'o2_abundance_20percent[22]',
                       'SN1a_2percent_density[23]', 'SN1a_4percent_density[24]','SN1a_6percent_density[25]',
                       'SN1a_8percent_density[26]', 'SN1a_10percent_density[27]','SN1a_15percent_density[28]',
                       'SN1a_20percent_density[29]']
        
        file_name = '%s/o2abundance%s.ascii' % (out_dir, current_scale)
        data = Table(ascii_data, names=ascii_names)
        ascii.write(data, output=file_name, overwrite=True) 

