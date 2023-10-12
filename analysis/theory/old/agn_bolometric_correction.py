


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


grid_dir = '/Users/sw376/Dropbox/Research/projects/synthesizer/tests/test_grid/'
grid_name = 'test_grid_blackholes'


grid = Grid(grid_name=grid_name, grid_dir=grid_dir)

spectra_type = 'total'

cmap = cmr.bubblegum
norm = mpl.colors.Normalize(vmin=4., vmax=7.) 

# initialise plot
fig = plt.figure(figsize=(3.5, 3.5))

left = 0.15
height = 0.8
bottom = 0.1
width = 0.8

# define main ax
ax = fig.add_axes((left, bottom, width, height))


# loop over log10ages

bolometric_correction = []

for i, log10T in enumerate(grid.log10T):
    
    grid_point = grid.get_grid_point((log10T,))
    sed = grid.get_spectra(grid_point)
    Lbol = sed.measure_bolometric_luminosity()
    Luv = np.interp(1500., sed.lam, sed.lnu) * (c / (1500 * Angstrom)).to('Hz').value
    bolometric_correction.append(Luv / Lbol)

ax.plot(grid.log10T, bolometric_correction, c='k', lw=2, label=r'$\rm FUV$')

# # inteprolation
# interp = interp1d(grid.log10T, bolometric_correction, kind='cubic')
# x = np.arange(4., 6, 0.01)
# ax.plot(x, interp(x), c='k', lw=2, label=r'$\rm FUV$')


ax.legend(fontsize=8, loc='upper left')
ax.set_xlabel(r'$\rm log_{10}(T/K)$')
ax.set_ylabel(r'$\rm L/L_{bol}$')


def kappa(Lbol, p):
    c1, k1, c2, k2 = p
    return c1*(Lbol/1E10)**-k1 + c1*(Lbol/1E10)**k2

for Lbol, lw in zip([10,11,12], [1,2,3]):

    # shen 2020 UV
    bc = kappa(10**Lbol, (1.862, -0.361, 4.870, -0.0063)) 
    ax.axhline(1/bc, c='k', lw=lw, alpha = 0.2)

    conversion = (1. * Lsun).to('erg/s')
    print(conversion)


    ax.text(4.05, 1/bc + 0.005, rf'$\rm Shen+2020\ [L_{{bol}}=10^{{{Lbol}}}\ L_{{\odot}}] \approx 3.8\times 10^{{{Lbol+33}}}\ erg/s$', fontsize=7, c='0.5')


ax.set_xlim([4.01, 5.99])
ax.set_ylim([0., 0.34])

filename = f'figs/agn_bolometric_correction.pdf'
print(filename)

fig.savefig(filename)
