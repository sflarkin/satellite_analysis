import yt
from yt.units import Mpc, Msun
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

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
            
            plt.savefig('halocenters/VELA07/%sid%s' % (VELA_a, halo[0]), bbox_inches='tight')
            #plt.show()
            
            plt.close()