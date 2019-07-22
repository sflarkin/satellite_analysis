import yt
from yt.units import Mpc, Msun
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def plotting_xyz_projection(ds, zoom, center, rvir, tomer_center, tomer_rvir, VELA_a, input_dir):
    #setup the figure grid to put the images onto
    fig = plt.figure()
    
    grid = AxesGrid(fig, (0.075,0.075,15,5),
    nrows_ncols = (6, 3),
    axes_pad = 0.05,
    label_mode = "L",
    share_all = True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%")
    
    ad = ds.all_data()
    par_type = 'darkmatter'
    masses = yt.np.unique(ad[(par_type, 'particle_mass')])
    for index in range(len(masses)):
        global filter_name, index1
        index1 = index
        filter_name = par_type + str(index)
        def mass_filter(pfilter, data):
            filter = data[(pfilter.filtered_type, 'particle_mass')] == masses[index]
            return filter
        yt.add_particle_filter(filter_name, function=mass_filter, filtered_type=par_type, requires=['particle_mass'])
        ds.add_particle_filter(filter_name)                       
                            
        global a,b,c
        a = yt.ParticlePlot(ds, (filter_name, 'particle_position_x'), (filter_name, 'particle_position_y'),\
                        (filter_name,'particle_mass'))
        a.set_unit((filter_name,'particle_mass'), 'Msun')
        a.set_figure_size(5)
        a.zoom(zoom)

        b = yt.ParticlePlot(ds, (filter_name, 'particle_position_y'), (filter_name, 'particle_position_z'),\
                        (filter_name,'particle_mass'))
        b.set_unit((filter_name,'particle_mass'), 'Msun')
        b.set_figure_size(5)
        b.zoom(zoom)

        c = yt.ParticlePlot(ds, (filter_name, 'particle_position_z'), (filter_name, 'particle_position_x'),\
                            (filter_name,'particle_mass'))
        c.set_unit((filter_name,'particle_mass'), 'Msun')
        c.set_figure_size(5)
        c.zoom(zoom)
    
        a.annotate_sphere(tomer_center, radius=(tomer_rvir, 'kpc'), circle_args={'color':'red'})
        a.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'black'})
        b.annotate_sphere(tomer_center, radius=(tomer_rvir, 'kpc'), circle_args={'color':'red'})
        b.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'black'})
        c.annotate_sphere(tomer_center, radius=(tomer_rvir, 'kpc'), circle_args={'color':'red'})
        c.annotate_sphere(center, radius=(rvir, 'kpc'), circle_args={'color':'black'})

        plot = a.plots[(filter_name, 'particle_mass')]
        plot.figure = fig
        plot.axes = grid[index*3].axes
        plot.cax = grid.cbar_axes[0]
        a._setup_plots()

        plot = b.plots[(filter_name, 'particle_mass')]
        plot.figure = fig
        plot.axes = grid[index*3+1].axes
        plot.cax = grid.cbar_axes[0]
        b._setup_plots()

        plot = c.plots[(filter_name, 'particle_mass')]
        plot.figure = fig
        plot.axes = grid[index*3+2].axes
        plot.cax = grid.cbar_axes[0]
        c._setup_plots()

    plt.savefig('%s/%sdarkmattercomparison.png' % (input_dir, VELA_a), bbox_inches='tight')
    plt.close()

def rvirprojections(ds, center, rvir, sp, out_dir, VELA_a, Id):
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
                          ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=sp)
    b = yt.ParticlePlot(ds, ('stars', 'particle_position_y'), ('stars', 'particle_position_z'),\
                          ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=sp)
    c = yt.ParticlePlot(ds, ('stars', 'particle_position_z'), ('stars', 'particle_position_x'),\
                          ('stars', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=sp)
            
    d = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_x'), ('darkmatter', 'particle_position_y'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=sp)
    e = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_y'), ('darkmatter', 'particle_position_z'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=sp)
    f = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_z'), ('darkmatter', 'particle_position_x'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*4, "kpc"), data_source=sp)
            
    a.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
    b.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
    c.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
    d.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
    e.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
    f.annotate_marker(center, 'x', plot_args={'s':100, 'color':'red'})
            
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
            
    plt.savefig('%s/%sid%s.png' % (out_dir, VELA_a, Id), bbox_inches='tight')
    plt.close()

    
    

def rvirprojectionsallhalos(ds, center, rvir, sp, tomer_center, tomer_rvir, halo_data, domain_width, input_dir, VELA_a):
    fig = plt.figure()
    
    grid = AxesGrid(fig, (0.075,0.075,15,5),
    nrows_ncols = (3, 3),
    axes_pad = 0.05,
    label_mode = "L",
    share_all = True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%")
    
    #create the 3 different x-y-z projection plots to add to the grid, one for stars and two for darkmatter
    a = yt.ParticlePlot(ds, ('stars', 'particle_position_x'), ('stars', 'particle_position_y'),\
                          ('stars', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    a.set_unit(('stars','particle_mass'), 'Msun')
    b = yt.ParticlePlot(ds, ('stars', 'particle_position_y'), ('stars', 'particle_position_z'),\
                          ('stars', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    b.set_unit(('stars','particle_mass'), 'Msun')
    c = yt.ParticlePlot(ds, ('stars', 'particle_position_z'), ('stars', 'particle_position_x'),\
                          ('stars', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    c.set_unit(('stars','particle_mass'), 'Msun')
            
    d = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_x'), ('darkmatter', 'particle_position_y'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    d.set_unit(('darkmatter','particle_mass'), 'Msun')
    e = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_y'), ('darkmatter', 'particle_position_z'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    e.set_unit(('darkmatter','particle_mass'), 'Msun')
    f = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_z'), ('darkmatter', 'particle_position_x'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    f.set_unit(('darkmatter','particle_mass'), 'Msun')
            
    g = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_x'), ('darkmatter', 'particle_position_y'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    g.set_unit(('darkmatter','particle_mass'), 'Msun')
    h = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_y'), ('darkmatter', 'particle_position_z'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    h.set_unit(('darkmatter','particle_mass'), 'Msun')
    i = yt.ParticlePlot(ds, ('darkmatter', 'particle_position_z'), ('darkmatter', 'particle_position_x'),\
                          ('darkmatter', 'particle_mass'), center=center, width=(rvir*3, "kpc/h"), data_source=sp)
    i.set_unit(('darkmatter','particle_mass'), 'Msun')
    
    #add the x marker for the largest halo from the rockstar catalog
    a.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    b.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    c.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    d.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    e.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    f.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    g.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    h.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    i.annotate_marker(center, 'o', plot_args={'s':100, 'color':'red'})
    
    #add the x marker for the tomer catalog
    a.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    b.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    c.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    d.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    e.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    f.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    g.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    h.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    i.annotate_marker(tomer_center, 'x', plot_args={'s':100, 'color':'black'})
    
    #add the circle for the tomer catalog for the 3rd row
    g.annotate_sphere(tomer_center, radius=(tomer_rvir, 'kpc/h'), circle_args={'color':'black'})
    h.annotate_sphere(tomer_center, radius=(tomer_rvir, 'kpc/h'), circle_args={'color':'black'})
    i.annotate_sphere(tomer_center, radius=(tomer_rvir, 'kpc/h'), circle_args={'color':'black'})
    
    #add the title for each column
    a.annotate_title('Column 1: X-Y Projection')
    b.annotate_title('Column 2: Y-Z Projection')
    c.annotate_title('Column 3: Z-X Projection')
    
    #add annotations for redshift and scale
    a.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, {'color':'black'})
    a.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    b.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, {'color':'black'})
    b.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    c.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, {'color':'black'})
    c.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    d.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, {'color':'black'})
    d.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    e.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, {'color':'black'})
    e.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    f.annotate_timestamp(corner='lower_left', redshift=True, time=False, draw_inset_box=False, {'color':'black'})
    f.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    g.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    h.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    i.annotate_scale(corner='lower_right', unit='Kpc/h', {'color':'black'})
    
    #add the annotations for what the dot and x is, and the color circles
    a.annotate_text([1,1], 'Red dot: Rockstar Center \n Black X: Nir Center', coord_system='axis'
    
    for halo in halo_data:
        x = float(halo[17])/domain_width
        y = float(halo[18])/domain_width
        z = float(halo[19])/domain_width
        rvir_smaller_halos = float(halo[11])*(float(VELA_a)/1000)
        center = [x, y, z]
        
        g.annotate_sphere(center, radius=(rvir_smaller_halos, 'kpc/h'), circle_args={'color':'red'})
        h.annotate_sphere(center, radius=(rvir_smaller_halos, 'kpc/h'), circle_args={'color':'red'})
        i.annotate_sphere(center, radius=(rvir_smaller_halos, 'kpc/h'), circle_args={'color':'red'})
            
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
        
    for letter in [g,h,i]:
        plot = letter.plots[('darkmatter', 'particle_mass')]
        plot.figure = fig
        plot.axes = grid[index].axes
        plot.cax = grid.cbar_axes[0]
        letter._setup_plots()
        index = index + 1
            
    plt.savefig('%s/%s.png' % (input_dir, VELA_a), bbox_inches='tight')
    plt.close()