import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr
import h5py
import os
from unyt import c, Angstrom, Lsun
from synthesizer.grid import Grid
import flare.plt as fplt
from scipy.interpolate import interp1d


grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids/'
grid_name = 'bpass-2.2.1-bin_chabrier03-0.1,100.0_cloudy-c17.03'


grid = Grid(grid_name=grid_name, grid_dir=grid_dir)

spectra_type = 'total'

# cmap = cmr.bubblegum
# norm = mpl.colors.Normalize(vmin=4., vmax=7.) 

colours = cmr.take_cmap_colors('cmr.pepper', len(grid.metallicity))


# initialise plot
fig = plt.figure(figsize=(3.5, 3.5))

left = 0.15
height = 0.8
bottom = 0.1
width = 0.8

# define main ax
ax = fig.add_axes((left, bottom, width, height))


# loop over log10ages



for j, (Z, col) in enumerate(zip(grid.metallicity, colours)):

    uv = []
    lyc = []

    for i, log10age in enumerate(grid.log10age):
        sed = grid.get_spectra((i,j), spectra_id='total')
        Lbol = sed.measure_bolometric_luminosity()
        Luv = np.interp(1500., sed.lam, sed.lnu) * (c / (1500 * Angstrom)).to('Hz').value
        uv.append(Luv / Lbol)
        # lyc.append(grid.log10Q['HI'][i,j]/Lbol)
        # print(10**grid.l/Lbol)

    ax.plot(grid.log10age, uv, c=col, lw=1, label = rf'$\rm Z={Z} $')
    # ax.plot(grid.log10age, lyc, c=col, lw=1, ls='--')

    # inteprolation
    # interp = interp1d(grid.log10age, bolometric_correction, kind='linear')
    # x = np.arange(6., 10., 0.01)
    # ax.plot(x, interp(x), c=col, lw=1, label = rf'$\rm Z={Z} $')


ax.legend(fontsize=7, loc='upper right', labelspacing=0.0)
ax.set_xlabel(r'$\rm log_{10}(age/yr)$')
ax.set_ylabel(r'$\rm L/L_{bol}$')
ax.set_xlim([6.,10.])
ax.set_ylim([0., 0.99])



filename = f'figs/stellar_bolometric_correction.pdf'
print(filename)

fig.savefig(filename)
