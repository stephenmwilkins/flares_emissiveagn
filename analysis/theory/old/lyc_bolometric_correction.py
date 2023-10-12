


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
bottom = 0.1
width = 0.8

# define main ax
ax = fig.add_axes((left, bottom, width, height))
ax2 = ax.twiny()


# -----------------------------------------------------
# -----------------------------------------------------
# ------- AGN

grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids/'
grid_name = 'blackholes_cloudy_c17.03_log10TBB'

grid = Grid(grid_name=grid_name, grid_dir=grid_dir)


# loop over log10T
bolometric_correction = []
for i, log10T in enumerate(grid.log10T):
    
    grid_point = grid.get_grid_point((log10T,))
    sed = grid.get_spectra(grid_point, spectra_id='incident')
    Lbol = sed.measure_bolometric_luminosity()
    Lbol = 1 * erg / s
    bolometric_correction.append(10**grid.log10Q['HI'][i]/Lbol)

print(np.max(bolometric_correction))

ax.plot(grid.log10T, np.log10(bolometric_correction), c='k', lw=2, label=r'$\rm cloudy\ AGN$')


ax.legend(fontsize=8, loc='upper left')
ax.set_xlim([4.01, 5.99])
ax.set_ylim([8., 10.49])
ax.set_xlabel(r'$\rm log_{10}(T_{BB}/K)$')
ax.set_ylabel(r'$\rm log_{10}(\dot{N}_{LyC}/L_{bol})$')

# -----------------------------------------------------
# -----------------------------------------------------
# ------- Stars

grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids/'
grid_name = 'bpass-2.2.1-bin_chabrier03-0.1,100.0_cloudy-c17.03'

grid = Grid(grid_name=grid_name, grid_dir=grid_dir)

# colours for each metallicity
colours = cmr.take_cmap_colors('cmr.amber', len(grid.metallicity))

for j, (Z, col) in enumerate(zip(grid.metallicity, colours)):

    lyc = []

    for i, log10age in enumerate(grid.log10age):
        sed = grid.get_spectra((i,j), spectra_id='incident')
        Lbol = sed.measure_bolometric_luminosity()
        # Lbol = np.interp(1500., sed.lam, sed.lnu) 
        lyc.append(10**grid.log10Q['HI'][i,j]/Lbol)

    print(min(np.log10(lyc)), max(np.log10(lyc)))

    ax2.plot(grid.log10age, np.log10(lyc), c=col, lw=1, label = rf'$\rm Z={Z} $')


ax2.legend(fontsize=7, loc='lower right', labelspacing=0.0)
ax2.set_xlim([6.01, 9.99])
ax2.set_xlabel(r'$\rm log_{10}(age/yr)$')


# -----------------------------------------------------
# -----------------------------------------------------

filename = f'figs/lyc_bolometric_correction.pdf'
print(filename)

fig.savefig(filename)
