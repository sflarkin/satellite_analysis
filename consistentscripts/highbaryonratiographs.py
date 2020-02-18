import yt
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
from mpl_toolkits.axes_grid1 import AxesGrid
from satellite_analysis.catalogreaders import consistentcatalogreader as consistent
from satellite_analysis.catalogreaders import tomercatalogreader as tomer
from satellite_analysis.graphs import stellarmassrelationfindingrvirprojections as rvirprj


#need to add arguments for VELA_dir and input_dir for the rockstarcatalogreader



input_dir = '/nobackupp2/sflarkin/GEN6/rockstar001analysis/VELA{}/hlists'
VELA_dir = '/nobackupp2/sflarkin/GEN6/VELA'
remove_subhalos = 'True'
out_dir = '/nobackupp2/sflarkin/GEN6/highbaryon'

halos_of_interest = [(500, '08', 758123.0)]#, (500, '08', 755273.0), (500, '08', 755429.0), (500, '08', 757482.0), (500, '08', 758089.0), (500, '08', 758039.0), (500, '08', 758194.0), (500, '08', 757452.0), (500, '08', 755133.0), (500, '08', 754857.0), (500, '08', 755204.0), (500, '10', 646850.0), (500, '10', 646757.0), (500, '10', 646891.0), (500, '10', 646952.0), (500, '10', 646730.0), (500, '11', 427964.0), (500, '11', 427983.0), (500, '14', 346146.0), (500, '14', 345761.0), (500, '14', 345762.0), (500, '15', 237964.0)]

VELAS_of_interest = np.unique([x[1] for x in halos_of_interest])

for VELA_numbers in VELAS_of_interest:
    VELA_dir = '/nobackupp2/sflarkin/GEN6/VELA{}'.format(VELA_numbers)
    input_dir = '/nobackupp2/sflarkin/GEN6/rockstar001analysis/VELA{}/hlists'.format(VELA_numbers)
    print(VELA_dir)
    print(input_dir)
    print('')
    
    consistent.consistent_catalog_reader(input_dir, remove_subhalos=remove_subhalos)

    VELA_snaps = glob.glob(VELA_dir + '/10MpcBox*')

    #now find the scales of the halos of interest
    
    scales_of_interest = np.unique([x[0] for x in halos_of_interest])
    
    #now find the consistent indices and VELA_index that match the scales of interest
    
    for scale in scales_of_interest:
        VELA_snap = glob.glob(VELA_dir + '/10MpcBox*{}*'.format(scale))
        if VELA_snap == []:
            print('Could not find corresponding VELA 10Mpc File for snapshot {} in VELA {}'.format(scale, VELA_numbers))
        else:
            #find the index of the consistent_catalog that matches the wanted timestep
            ds = yt.load(VELA_snap[0])
            domain_width = float(ds.domain_width.in_units('Mpccm/h')[0])
            
            catalog_index = [x for x,i in enumerate(consistent.consistent_file_index) if i == str(scale)]
            print(catalog_index)
            
            halo_data = consistent.halo_data_sorted[catalog_index[0]]
            ids_of_interest = [int(x[2]) for x in halos_of_interest if x[1] == str(VELA_numbers)]
            
            for halos in halo_data:
                if int(halos[1]) in ids_of_interest:
                    x = float(halos[17])/domain_width
                    y = float(halos[18])/domain_width
                    z = float(halos[19])/domain_width
                    rvir = float(halos[11])*(float(scale)/1000) / .7
                    center = [x, y, z]
                    print(center, rvir)
                    sp = ds.sphere(center, (rvir*2, 'kpc'))
                    
                    graph_size = rvir*4
                    
                    #now plot the graph
                    
                    fig = plt.figure()
    
                    grid = AxesGrid(fig, (0.075,0.075,10,5),
                    nrows_ncols = (3, 2),
                    axes_pad = 1.0,
                    label_mode = "L",
                    share_all = False,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="3%",
                    cbar_pad="0%")
    
                    #create the 3 different x-y-z projection plots to add to the grid, one for stars and two for darkmatter
                    a = yt.ParticlePlot(ds, ('stars', 'particle_position_x'), ('stars', 'particle_position_y'),\
                                          ('stars', 'particle_mass'), center=center, width=(graph_size, "kpc"), data_source=sp)
                    a.set_unit(('stars','particle_mass'), 'Msun')
                    b = yt.ParticlePlot(ds, ('stars', 'particle_position_y'), ('stars', 'particle_position_z'),\
                                          ('stars', 'particle_mass'), center=center, width=(graph_size, "kpc"), data_source=sp)
                    b.set_unit(('stars','particle_mass'), 'Msun')
                    c = yt.ParticlePlot(ds, ('stars', 'particle_position_z'), ('stars', 'particle_position_x'),\
                                          ('stars', 'particle_mass'), center=center, width=(graph_size, "kpc"), data_source=sp)
                    c.set_unit(('stars','particle_mass'), 'Msun')
            
            
                    g = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_x'), ('darkmatter', 'particle_position_y'),\
                                          ('darkmatter', 'particle_mass'), center=center, width=(graph_size, "kpc"), data_source=sp)
                    g.set_unit(('darkmatter','particle_mass'), 'Msun')
                    h = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_y'), ('darkmatter', 'particle_position_z'),\
                                          ('darkmatter', 'particle_mass'), center=center, width=(graph_size, "kpc"), data_source=sp)
                    h.set_unit(('darkmatter','particle_mass'), 'Msun')
                    i = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_z'), ('darkmatter', 'particle_position_x'),\
                                          ('darkmatter', 'particle_mass'), center=center, width=(graph_size, "kpc"), data_source=sp)
                    i.set_unit(('darkmatter','particle_mass'), 'Msun')
    
                    #add the x marker for the largest halo from the rockstar catalog
                    a.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                    b.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                    c.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                    
                    g.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                    h.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                    i.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
                    
                    #add the circle
                    a.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'red'})
                    b.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'red'})
                    c.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'red'})
                            
                    g.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'red'})
                    h.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'red'})
                    i.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'red'})
    
                    #add the titles 
                    a.annotate_title('Star Particles: X-Y Projection')
                    b.annotate_title('Star Particles: Y-Z Projection')
                    c.annotate_title('Star Particles: Z-X Projection')
                    
                    g.annotate_title('Smallest Dark Matter Particles: X-Y Projection')
                    h.annotate_title('Smallest Dark Matter Particles: Y-Z Projection')
                    i.annotate_title('Smallest Dark Matter Particles: Z-X Projection')
    
                    #add annotations for redshift and scale
                    a.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, text_args={'color':'black'})
                    a.annotate_scale(corner='lower_right', unit='kpc', coeff=10, size_bar_args={'color': 'black'})
                    b.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, text_args={'color':'black'})
                    b.annotate_scale(corner='lower_right', unit='kpc', coeff=10, size_bar_args={'color': 'black'})
                    c.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, text_args={'color':'black'})
                    c.annotate_scale(corner='lower_right', unit='kpc', coeff=10, size_bar_args={'color': 'black'})
                    
                    g.annotate_scale(corner='lower_right', unit='kpc', coeff=10, size_bar_args={'color': 'black'})
                    h.annotate_scale(corner='lower_right', unit='kpc', coeff=10, size_bar_args={'color': 'black'})
                    i.annotate_scale(corner='lower_right', unit='kpc', coeff=10, size_bar_args={'color': 'black'})
    
                    #now find the halos close to this halo to plot circles over
                    for close_halos in halo_data:
                        x = float(close_halos[17])/domain_width
                        y = float(close_halos[18])/domain_width
                        z = float(close_halos[19])/domain_width
                        rvir_smaller_halos = float(close_halos[11])*(float(scale)/1000) / .7
                        center_close = [x, y, z]
                        distance_to_center = ((center[0]*center_close[0])**2 +  (center[1]*center_close[1])**2 + (center[2]*center_close[2])**2)**(1/2) * domain_width
                        
                        #if the halo is within the graph, add its circle to the grids
                        if distance_to_center < rvir*2 and halos[1] != close_halos[1]:
                            a.annotate_sphere(center_close, radius=(rvir_smaller_halos, 'kpc'), circle_args={'color':'red'})
                            b.annotate_sphere(center_close, radius=(rvir_smaller_halos, 'kpc'), circle_args={'color':'red'})
                            c.annotate_sphere(center_close, radius=(rvir_smaller_halos, 'kpc'), circle_args={'color':'red'})
                            
                            g.annotate_sphere(center_close, radius=(rvir_smaller_halos, 'kpc'), circle_args={'color':'red'})
                            h.annotate_sphere(center_close, radius=(rvir_smaller_halos, 'kpc'), circle_args={'color':'red'})
                            i.annotate_sphere(center_close, radius=(rvir_smaller_halos, 'kpc'), circle_args={'color':'red'})
            
                    index = 0
                    for letter in [a,b,c]:
                        plot = letter.plots[('stars', 'particle_mass')]
                        plot.figure = fig
                        plot.axes = grid[index].axes
                        plot.cax = grid.cbar_axes[index]
                        letter._setup_plots()
                        index = index + 2
                    
                    index = 1
                    for letter in [g,h,i]:
                        plot = letter.plots[('darkmatter', 'particle_mass')]
                        plot.figure = fig
                        plot.axes = grid[index].axes
                        plot.cax = grid.cbar_axes[index]
                        letter._setup_plots()
                        index = index + 2
            
                    plt.savefig('{}/VELA{}_Scale{}_id{}.png'.format(out_dir, VELA_numbers, scale, halos[1]), bbox_inches='tight')
                    plt.close()   
            
            
            
            
            
            
                              