
import numpy as np
import matplotlib.cm as cm
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.colors import Normalize

import scipy.stats as stats
import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.limits
import flares_utility.plt
import flares_utility.analyse as analyse
import flares_utility.stats

from unyt import c, Msun, yr, Lsun

filename = '/Users/sw376/Dropbox/Research/data/simulations/flares/flares_no_particlesed.hdf5'

flares = analyse.analyse(filename, default_tags=False)

# flares.list_datasets()


tag = '010_z005p000' # z=5


quantities = []

quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mass',
                  'name': 'BH_Mass', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mdot',
                  'name': 'BH_Mdot', 'log10': True})

D = {}
s = {}

# --- get quantities (and weights and deltas)
D = flares.get_datasets(tag, quantities)


xlimits = Mbh_limit = [6.51, 9.49] # Msun
ylimits = Mbhdot_limit = np.array([-5.99, 0.99]) # Msun/yr

conversion = 0.1 * Msun * c**2 / yr
log10conversion = np.log10(conversion.to('erg/s').value)

log10Lbol = D['log10BH_Mdot'] + log10conversion


fig = plt.figure(figsize = (3.5, 2.5))

left  = 0.15
height = 0.7
bottom = 0.2
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

s = D['log10BH_Mass']>6.5


log10T = np.log10(2.24E9 * D['BH_Mdot'][s] ** (1 / 4) * D['BH_Mass'][s]**-0.5)

# all BHs >10^8 Msol
N, bin_edges = np.histogram(log10T, bins = np.arange(4., 6.5, 0.1))
# ax.step(bin_edges[:-1], N, where='post', c='k', label = r'$\rm M_{\bullet}>10^{7}$', )
ax.stairs(N, bin_edges, color='0.9', label = r'$\rm M_{\bullet}>10^{7}$', fill=True )


# only BHs >10^8 Msol AND Lbol>1E45
s2 = log10Lbol[s] > 45.

N, bin_edges = np.histogram(log10T[s2], bins = np.arange(4., 6.5, 0.1))
ax.step(bin_edges[:-1], N, where='post', c='k', label = r'$\rm M_{\bullet}>10^{7}; L_{bol}>10^{45}\ erg/s$')



ax.set_xlim([4.01, 5.99])
# ax.set_ylim(Mbhdot_limit)

ax.legend(fontsize=8, loc='upper left')
ax.set_xlabel(r'$\rm log_{10}(T/K)$')
ax.set_ylabel(r'$\rm N$')


filename = f'figs/Tdist.pdf'
print(filename)

fig.savefig(filename)

