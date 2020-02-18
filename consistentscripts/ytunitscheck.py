from yt.utilities.cosmology import Cosmology
co = Cosmology(hubble_constant = .70)
from yt.units.yt_array import YTQuantity
import numpy as np

z = 2
a = 1/3

critical_density = co.critical_density(z)
print('critical_density: {}'.format(critical_density))
print("critical_density.in_units('Msun*h**2/kpc**3'): {}".format(critical_density.in_units('Msun*h**2/kpc**3')))

#units we want to convert
Mvir = YTQuantity(1279000000000.0, 'Msun/h')
Rvir_comoving = YTQuantity(279.501, 'kpc/h')
#these Rvirs are in comoving distance, so we multiple by a to convert to proper
Rvir_proper = Rvir_comoving*a

volume = (4/3)*np.pi*(Rvir_proper**3)
density = Mvir/volume
print('density: {}'.format(density))
print('')


print('Normal division, no unit conversion: {}'.format(density/critical_density))
print('Normal division, unit conversion to be same as density: {}'.format(density/(critical_density.in_units('Msun*h**2/kpc**3'))))
print('')
print('Now I compare this to the division of just their floats, or with the units removed.')
print('Division of their floats: {}'.format(float(density)/float(critical_density.in_units('Msun*h**2/kpc**3'))))
print('')
print('There is a factor of .7**2 diference between them, shown below.')
print(density/critical_density*(.7**2))
print('This seems wrong to me, as the units divided are the same, so they should cancel out, yet the h**2 seems to remain.')


