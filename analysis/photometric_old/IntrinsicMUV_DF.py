
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.lines import Line2D

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.analyse as analyse

from synthesizer.utils import Lnu_to_M

from unyt import c, Msun, yr, g, s, Angstrom, Hz

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
redshifts = [10, 9, 8, 7, 6, 5]


binw = 1.
X_limits = [-18.01, -25.99]
Y_limits = [-7.9, -2.01]


filename = '/Users/sw376/Dropbox/Research/data/simulations/flares/flares_no_particlesed.hdf5'
flares = analyse.analyse(filename)

obs_colours = cmr.take_cmap_colors('cmr.infinity', 5, cmap_range = (0.1,0.9))



N = len(redshifts)
left = 0.1
top = 0.95
bottom = 0.15
right = 0.975
panel_width = (right-left)/N
panel_height = top-bottom
fig, axes = plt.subplots(2, 3, figsize = (7,5), sharey = True, sharex = True)
plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.1)




# --- FLARES


# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3


bin_edges = np.arange(*[-27.,-16.], binw)
bin_centres = bin_edges[:-1]+binw/2


# converting MBHacc units to M_sol/yr
h = 0.6777  # Hubble parameter
bhacc_conversion = h * 6.446E23  * g / s
conversion = 0.1 * bhacc_conversion * c**2 # Lbol
# conversion /= 9.0 # convert to LUV
# conversion /= (c/(1500. * Angstrom))
conversion *= 7.9E-17 / Hz
conversion = conversion.to('erg/s/Hz').value
log10conversion = np.log10(conversion)

print(log10conversion)


for z, c, ax in zip(redshifts, cmr.take_cmap_colors('cmr.bubblegum_r', len(redshifts)), axes.flatten()):

    tag = flares.tag_from_zed[z]

    # get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))


    # Lbol BH

    BH_Mdot = flares.load_dataset(tag, *['Galaxy', 'BH_Mdot'])

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.array(BH_Mdot[sim]) * conversion
        x = Lnu_to_M(x)

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w

    # ax.plot(bin_centres[N>4], np.log10(phi[N>4]), ls = ':', c=c, lw=1)
    ax.plot(bin_centres, np.log10(phi), ls = '-', c=c, lw=2, alpha=0.3)
    ax.plot(bin_centres[N>4], np.log10(phi[N>4]), ls = '-', c=c, lw=2)



    # Lbol * 

    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))

    LUVStar = flares.load_dataset(tag, *['Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Intrinsic', 'FUV'])

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.array(LUVStar[sim])
        x = x[x>0.0]
        x = Lnu_to_M(x)

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w

    # ax.plot(bin_centres[N>4], np.log10(phi[N>4]), ls = '--', c=c, lw=1)
    ax.plot(bin_centres, np.log10(phi), ls = '--', c=c, lw=1, alpha=0.3)
    ax.plot(bin_centres[N>4], np.log10(phi[N>4]), ls = '--', c=c, lw=1)



    # TOTAL

    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.array(LUVStar[sim]) + np.array(BH_Mdot[sim]) * conversion
        x = x[x>0.0]
        x = Lnu_to_M(x)

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w


    # ax.plot(bin_centres[N>4], np.log10(phi[N>4]), ls = '-', c=c, lw=2, alpha=0.3)
    ax.plot(bin_centres, np.log10(phi), ls = ':', c=c, lw=1, alpha=0.3)
    ax.plot(bin_centres[N>4], np.log10(phi[N>4]), ls = ':', c=c, lw=1)




    ax.text(0.5, 1.02, rf'$\rm z={z}$', horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize = 7)

    ax.set_xlim(X_limits)
    ax.set_ylim(Y_limits)

    # observations

    # Maiolino
    if z==5:

        # c = obs_colours[0]
        c = 'k'
        x = [-18., -19.5]
        y = [-3.06, -3.83]
        xerr = [[0.75, 0.75],[0.75, 0.75]]
        yerr = [[0.26,0.44],[0.18,0.23]]

        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o", markersize=4, lw=1, color=c, label=r'$\rm Maiolino+23$')


    # Kocevski
    if z==5:

        # c = obs_colours[1]
        c = 'k'
        x = [-19.4]
        y = [-4.97]
        yerr = [[0.3],[0.2]]

        ax.errorbar(x, y, yerr=yerr, fmt="D", markersize=4, lw=1, color=c, label=r'$\rm Kocevski+23$')

    # Matthe
    if z==5:

        # c = obs_colours[2]
        c = 'k'
        x = [-20, -19, -18]
        y = [-5.08, -4.8, -4.98]
        xerr = [0.5, 0.5, 0.5]
        yerr = [[0.26, 0.18, 0.23],[0.16, 0.13, 0.15]]

        ax.errorbar(x, y, yerr=yerr, fmt="s", markersize=4, lw=1, color=c, label=r'$\rm Matthee+23$')

    ax.legend(fontsize=7, labelspacing=0.2)

fig.text(0.04, bottom+(top-bottom)/2, r'$\rm\log_{10}(\phi/Mpc^{-3}\ dex^{-1})$', rotation = 90, va='center', fontsize = 9)
fig.text(left+(right-left)/2, 0.08, r'$\rm M_{1500}$', ha='center', fontsize = 9)


handles = []
handles.append(Line2D([0], [0], label=r'$\rm stars$', color='0.5', lw=1, alpha=0.5, ls='--', c='k'))
handles.append(Line2D([0], [0], label=r'$\rm blackholes$', color='0.5', lw=2, alpha=0.5, ls='-', c='k'))
handles.append(Line2D([0], [0], label=r'$\rm total$', color='0.5', lw=1, alpha=0.5, ls=':', c='k'))

fig.legend(handles=handles, fontsize=8, labelspacing=0.1, loc = 'outside lower center', ncol=len(handles))

fig.savefig(f'figs/IntrinsicMUV_DF.pdf')


fig.clf()
