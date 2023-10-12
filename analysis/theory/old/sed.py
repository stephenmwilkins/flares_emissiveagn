


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr
import h5py
import os
from unyt import c, Angstrom, Lsun, erg, s
from synthesizer.grid import Grid
import flare.plt as fplt
from scipy.interpolate import interp1d






# initialise plot
fig = plt.figure(figsize=(3.5, 3.5))

left = 0.15
height = 0.8
bottom = 0.15
width = 0.8

# define main ax
ax = fig.add_axes((left, bottom, width, height))




grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids/'
grid_name = 'blackholes_cloudy_c17.03_log10TBB'

grid = Grid(grid_name=grid_name, grid_dir=grid_dir)

log10T = 5.5
    
grid_point = grid.get_grid_point((log10T,))

for spectra_id, lw, alpha in zip(['total', 'incident'], [2,1], [0.4, 1.0]):

    sed = grid.get_spectra(grid_point, spectra_id=spectra_id)
    Lbol = sed.measure_bolometric_luminosity(method='quad')
    print(Lbol)
    ax.plot(np.log10(sed.lam)-4., np.log10(sed.nu*sed.lnu)+45, c='k', lw=lw, alpha=alpha, label = rf'$\rm {spectra_id}$')

    Luv = sed.measure_window_lnu([1400,1600]) * (c / (1500 * Angstrom)).to('Hz').value
    print(np.log10(Lbol), np.log10(Luv))


ax.axvline(np.log10(0.15), c='k', alpha=0.1, lw=1)

ax.legend(fontsize=8, loc='upper left')
ax.set_xlim([-4., 4.])
ax.set_xlabel(r'$\rm log_{10}(\lambda/\mu m)$')

ax.set_ylim([40., 45.99])
ax.set_ylabel(r'$\rm log_{10}(\nu L_{\nu}/erg\ s^{-1})$')


# ax.set_ylim([28., 33.])
# ax.set_ylabel(r'$\rm log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1})$')

filename = f'figs/sed.pdf'
print(filename)

fig.savefig(filename)
