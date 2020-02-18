#this is for the snapnumber3x2

import yt
import glob
import numpy as np
import random
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
plt.switch_backend('agg')
from satellite_analysis.catalogreaders import consistentcatalogreader as consistent 
from satellite_analysis.catalogreaders import tomercatalogreader as tomer

def parse():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir')
    parser.add_argument('VELA_dir')
    parser.add_argument('VELA_number')
    parser.add_argument('out_dir')
    args = vars(parser.parse_args())
    return args

args = parse()

input_dir = args['input_dir']
VELA_dir = args['VELA_dir']
VELA_number = args['VELA_number']
out_dir = args['out_dir']

#load in the rockstar catatalog halo data
consistent.consistent_catalog_reader(input_dir, halomass=1e+07)

#find the VELA snapshots
VELA_snaps = glob.glob(VELA_dir + '/10MpcBox_csf512_a0*')
VELA_index = []
for snap in VELA_snaps:
    period = [pos for pos, char in enumerate(snap) if char == '.']
    number = snap[period[-2]+1:period[-1]]
    VELA_index.append(number)
VELA_snaps.sort()
VELA_index.sort()

for index in consistent.snapshot_index:
    VELA_a = consistent.consistent_file_index[index]
    print('Generating Graph Grid for snapshot:', VELA_a)
    position = [pos for pos, loc in enumerate(VELA_index) if loc == VELA_a]
    if position == [] or len(position) > 1:
        print('Could not find corresponding VELA 10Mpc File for snapshot:', VELA_a)
    else:
        ds = yt.load(VELA_snaps[position[0]])
        domain_width = float(ds.domain_width.in_units('Mpc/h')[0])
        ad = ds.all_data()
        masses = yt.np.unique(ad[('darkmatter', 'particle_mass')])
        #filter out the darkmatter0 particles
        def mass_filter(pfilter, data):
            filter = data[(pfilter.filtered_type, 'particle_mass')] == masses[0]
            return filter
        yt.add_particle_filter('darkmatter0', function=mass_filter, filtered_type='darkmatter', requires=['particle_mass'])
        ds.add_particle_filter('darkmatter0')
        
        scale = float(ds.scale_factor)

        fig = plt.figure()
    
        grid = AxesGrid(fig, (0.075,0.075,10,5),
        nrows_ncols = (3, 1),
        axes_pad = 1.0,
        label_mode = "L",
        share_all = False,
        cbar_location="right",
        cbar_mode="each",
        cbar_size="3%",
        cbar_pad="0%")
        
        zoom = 10
        
        x0 = float(consistent.halo_data_sorted[index][0][17]) * scale / domain_width
        y0 = float(consistent.halo_data_sorted[index][0][18]) * scale / domain_width
        z0 = float(consistent.halo_data_sorted[index][0][19]) * scale / domain_width
        center = [x0, y0, z0]
        
        a = yt.ParticlePlot(ds, ('darkmatter0', 'particle_position_x'), ('darkmatter0', 'particle_position_y'),\
                              ('darkmatter0', 'particle_mass'), center=center)
        a.set_unit(('darkmatter0','particle_mass'), 'Msun')
        a.zoom(zoom)
        b = yt.ParticlePlot(ds, ('darkmatter0', 'particle_position_y'), ('darkmatter0', 'particle_position_z'),\
                              ('darkmatter0', 'particle_mass'), center=center)
        b.set_unit(('darkmatter0','particle_mass'), 'Msun')
        b.zoom(zoom)
        c = yt.ParticlePlot(ds, ('darkmatter0', 'particle_position_z'), ('darkmatter0', 'particle_position_x'),\
                              ('darkmatter0', 'particle_mass'), center=center)
        c.set_unit(('darkmatter0','particle_mass'), 'Msun')
        c.zoom(zoom)
        

        for halos in consistent.halo_data_sorted[index]:
            x = float(halos[17]) * scale / domain_width
            y = float(halos[18]) * scale / domain_width
            z = float(halos[19]) * scale / domain_width
            r = float(halos[11]) * scale / .7
            center = [x, y, z]
            a.annotate_sphere(center, radius=(r, 'kpc'), circle_args={'color':'red'})
            b.annotate_sphere(center, radius=(r, 'kpc'), circle_args={'color':'red'})
            c.annotate_sphere(center, radius=(r, 'kpc'), circle_args={'color':'red'})
                            
    
        index = 0
        for letter in [a,b,c]:
            plot = letter.plots[('darkmatter0', 'particle_mass')]
            plot.figure = fig
            plot.axes = grid[index].axes
            plot.cax = grid.cbar_axes[index]
            letter._setup_plots()
            index = index + 1

        plt.savefig('{}/CatalogProjectionYT_VELA{}_Scale{}.png'.format(out_dir, VELA_number, str(scale)[2:5]), bbox_inches='tight')
        plt.close()   
    
        