
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.analyse as analyse
from flares_utility.plt import measure_weighted_median, weighted_median

from unyt import c, Msun, yr, g, s, Angstrom

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
redshifts = [10, 9, 8, 7, 6, 5]
xlimits = [54.51, 56.99]
ylimits = [-3.9, 1.99]


filename = '/Users/sw376/Dropbox/Research/data/simulations/flares/flares_no_particlesed.hdf5'
flares = analyse.analyse(filename, default_tags=False)
flares.list_datasets()

filename = '/Users/sw376/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v4.hdf5'
flares2 = analyse.analyse(filename, default_tags=False)
flares2.list_datasets()

N = len(redshifts)
left = 0.1
top = 0.95
bottom = 0.15
right = 0.9
panel_width = (right-left)/N
panel_height = top-bottom
fig, axes = plt.subplots(2, 3, figsize = (7,5), sharey = True, sharex = True)
plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.1)

cax = fig.add_axes([right, bottom, 0.03, top-bottom])


# converting MBHacc units to M_sol/yr

conversion = 0.1 * Msun * c**2 / yr
log10bolometric_conversion = np.log10(conversion.to('erg/s').value)
log10conversion = log10bolometric_conversion + 10. # approximate erg/s -> 1/s



quantities = []
quantities.append({'path': 'Galaxy', 'dataset': f'Mstar_30',
                  'name': 'Mstar', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mass',
                  'name': 'BH_Mass', 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': f'BH_Mdot',
                  'name': 'BH_Mdot', 'log10': True})
quantities.append({'path': 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Intrinsic/', 'dataset': f'FUV',
                  'name': 'FUV_Intrinsic', 'log10': True})

quantities2 = []
quantities2.append({'path': 'Galaxy', 'dataset': f'IonisingEmissivity',
                  'name': 'IonisingEmissivity', 'log10': True})




norm = Normalize(vmin=8., vmax=11.5)
cmap = cmr.lavender


for z, c, ax in zip(redshifts, cmr.take_cmap_colors('cmr.bubblegum_r', len(redshifts)), axes.flatten()):

    # ax.fill_between([0, 45.],[-9,-9],[3,3],color='k',alpha=0.05)
    ax.axhline(0.0, color='k', lw=2, alpha = 0.1)

    tag = flares.tag_from_zed[z]

    D = flares.get_datasets(tag, quantities)
    D2 = flares2.get_datasets(tag, quantities2)

    log10Mstar = D['log10Mstar']
    log10Lbhs = D['log10BH_Mdot'] + log10conversion # should be erg/s/Hz
    log10Lbhs_bol = D['log10BH_Mdot'] + log10bolometric_conversion # should be erg/s
    log10Lstars = D2['log10IonisingEmissivity']

    log10Ltotal = np.log10(10**log10Lbhs + 10**log10Lstars)

    ratio = log10Lbhs - log10Lstars

    # everything
    ax.scatter(log10Ltotal, ratio, s=1, c='k', zorder=1, alpha=0.2)

    s = log10Lbhs_bol > 45.
    s2 = s & (D['log10BH_Mass']>7.)

    ax.scatter(log10Ltotal[s2], ratio[s2], s=3, c='k', zorder=1)
    ax.scatter(log10Ltotal[s2], ratio[s2], s=1, c=cmap(norm(log10Mstar[s2])), zorder=2)

    ax.scatter(log10Ltotal[s], ratio[s], s=1, c=cmap(norm(log10Mstar[s])), zorder=2, alpha=0.5)

    


    bins = np.arange(44, 47., 0.25)
    # N, bins, bincen, cuts, med = measure_weighted_median(log10Ltotal[s], ratio[s], D['weight'][s], bins=bins, weighted=False)

    # ax.plot(bincen, med, c='k')
    
    weighted_median(ax, log10Ltotal[s], ratio[s], D['weight'][s], bins=bins, weighted=False)


    ax.text(0.5, 1.02, rf'$\rm z={z}$', horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize = 7)
    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)



fig.text(0.04, bottom+(top-bottom)/2, r'$\rm \log_{10}(\dot{N}_{\bullet,\ LyC}/\dot{N}_{\star,\ LyC})$', rotation = 90, va='center', fontsize = 9)
fig.text(left+(right-left)/2, 0.08, r'$\rm\log_{10}(\dot{N}_{\bullet,\ LyC}+\dot{N}_{\star,\ LyC}/s^{-1})$', ha='center', fontsize = 9)


# add colourbar
cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(cmapper, cax=cax, orientation='vertical')
cax.yaxis.tick_right()
cax.yaxis.set_label_position('right')
cax.set_ylabel(r'$\rm log_{10}(M_{\star}/M_{\odot})$', fontsize=8)
cax.tick_params(axis='y', labelsize=6)



fig.savefig(f'figs/relative_LyC.pdf')
fig.savefig(f'figs/relative_LyC.png')

fig.clf()