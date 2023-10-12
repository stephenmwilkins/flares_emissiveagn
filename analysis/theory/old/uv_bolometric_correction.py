


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



# spectra type to use for both models
spectra_id = 'total' 

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
    sed = grid.get_spectra(grid_point, spectra_id=spectra_id)
    # Lbol = sed.measure_bolometric_luminosity()
    Lbol = 1 * erg/s

    Luv = sed.measure_window_lnu([1400,1600]) * (c / (1500 * Angstrom)).to('Hz').value
    # Luv = np.interp(1500., sed.lam, sed.lnu) * (c / (1500 * Angstrom)).to('Hz').value
    bolometric_correction.append(Luv / Lbol)

ax.plot(grid.log10T, bolometric_correction, c='k', lw=2, label=r'$\rm cloudy\ AGN$')


# add Shen corrections
def kappa(Lbol, p):
    c1, k1, c2, k2 = p
    return c1*(Lbol/1E10)**-k1 + c1*(Lbol/1E10)**k2

for Lbol, lw in zip([10,11,12], [1,2,3]):

    bc = kappa(10**Lbol, (1.862, -0.361, 4.870, -0.0063)) 
    ax.axhline(1/bc, c='k', lw=lw, alpha = 0.2)

    conversion = (1. * Lsun).to('erg/s')
    ax.text(4.05, 1/bc + 0.01, rf'$\rm Shen+2020\ [L_{{bol}}=10^{{{Lbol}}}\ L_{{\odot}}] \approx 3.8\times 10^{{{Lbol+33}}}\ erg/s$', fontsize=7, c='0.5')


ax.legend(fontsize=8, loc='upper left')
ax.set_xlim([4.01, 5.99])
ax.set_ylim([0., 0.99])
ax.set_xlabel(r'$\rm log_{10}(T_{BB}/K)$')
ax.set_ylabel(r'$\rm L/L_{bol}$')

# -----------------------------------------------------
# -----------------------------------------------------
# ------- Stars

grid_dir = '/Users/sw376/Dropbox/Research/data/synthesizer/grids/'
grid_name = 'bpass-2.2.1-bin_chabrier03-0.1,100.0_cloudy-c17.03'

grid = Grid(grid_name=grid_name, grid_dir=grid_dir)

# colours for each metallicity
colours = cmr.take_cmap_colors('cmr.amber', len(grid.metallicity))

for j, (Z, col) in enumerate(zip(grid.metallicity, colours)):

    uv = []
    lyc = []

    for i, log10age in enumerate(grid.log10age):
        sed = grid.get_spectra((i,j), spectra_id=spectra_id)
        Lbol = sed.measure_bolometric_luminosity()
        Luv = np.interp(1500., sed.lam, sed.lnu) * (c / (1500 * Angstrom)).to('Hz').value
        uv.append(Luv / Lbol)
        # lyc.append(grid.log10Q['HI'][i,j]/Lbol)
        # print(10**grid.l/Lbol)

    ax2.plot(grid.log10age, uv, c=col, lw=1, label = rf'$\rm Z={Z} $')


ax2.legend(fontsize=7, loc='upper right', labelspacing=0.0)
ax2.set_xlim([6.01, 9.99])
ax2.set_xlabel(r'$\rm log_{10}(age/yr)$')


# -----------------------------------------------------
# -----------------------------------------------------

filename = f'figs/uv_bolometric_correction.pdf'
print(filename)

fig.savefig(filename)
