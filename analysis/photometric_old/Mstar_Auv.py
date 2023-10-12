
import numpy as np
import matplotlib.cm as cm
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.colors import Normalize
import scipy.stats as stats
from scipy.stats import pearsonr

import cmasher as cmr

import h5py


import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.limits
import flares_utility.plt
import flares_utility.analyse as analyse
import flares_utility.stats

filename = '/Users/sw376/Dropbox/Research/data/simulations/flares/flares_no_particlesed.hdf5'

flares = analyse.analyse(filename, default_tags=False)

# flares.list_datasets()



tag = '010_z005p000'


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': f'Mstar_30',
                  'name': 'Mstar', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mass',
                  'name': 'BH_Mass', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mdot',
                  'name': 'BH_Mdot', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_los',
                  'name': 'BH_los', 'log10': True})



fig = plt.figure(figsize = (3.5, 3.5))

left  = 0.15
height = 0.8
bottom = 0.15
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


x = 'Mstar'
y = 'BH_Mass'
z = 'BH_Mdot'


xlimits = np.array([9., 11.])
ylimits = np.array([0., 50.])



norm = Normalize(vmin=9., vmax=11.)
cmap = cm.plasma

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': f'Mstar_30',
                  'name': 'Mstar', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mass',
                  'name': 'BH_Mass', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_los',
                  'name': 'BH_los', 'log10': True})
quantities.append({'path': 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Intrinsic/', 'dataset': f'FUV',
                  'name': 'FUV_Intrinsic', 'log10': False})
quantities.append({'path': 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/', 'dataset': f'FUV',
                  'name': 'FUV_DustModelI', 'log10': False})
quantities.append({'path': 'Galaxy', 'dataset': f'DTM',
                  'name': 'DTM', 'log10': False})

D = flares.get_datasets(tag, quantities)



Astar = 2.5*np.log10(D['FUV_Intrinsic']/D['FUV_DustModelI'])

kappa = 0.0795
gamma = -1.0
Abh = kappa * D['DTM'] * D['BH_los'] * (1500./5500.)**gamma

r = pearsonr(Astar, Abh)
print(r)

s = D['log10Mstar']>9.

print(np.median(D['DTM'][s]))
print(np.median(D['BH_los'][s]))
print(np.median(Abh[s]))

ax.scatter(D['log10Mstar'][s], Abh[s], s=1, c='r', label = 'bhs', alpha = 0.2)
ax.scatter(D['log10Mstar'][s], Astar[s], s=1, c='b', label = 'stars', alpha = 0.2)

ax.set_ylim(ylimits)
ax.set_xlim(xlimits)

ax.set_xlabel(r'$\rm log_{10}(M_{\star}/M_{\odot})$')
ax.set_ylabel(r'$\rm A_{FUV}$')

ax.legend()

fig.savefig(f'figs/Mstar_Auv.pdf')
