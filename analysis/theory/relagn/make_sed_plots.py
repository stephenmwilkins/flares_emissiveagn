
import numpy as np
import matplotlib.pyplot as plt
from relagn import relagn
from unyt import c
import cmasher as cmr

# set style
plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')


variants = [
    ('M', [7,8,9], r'log_{10}(M_{\bullet}/M_{\odot})', 'cmr.torch', True),
    ('a', [0.0, 0.5, 0.9, 0.99, 0.998], r'a', 'cmr.cosmic', False),
    ('log_mdot', [-3, -2, -1, 0], r'log_{10}(\dot{M}/\dot{M}_{Edd})', 'cmr.ember', False),
    ('cos_inc', [0.09, 0.5, 0.98], r'\cos(\theta)', 'cmr.bubblegum', False),
]



for variant in variants: 

    fig = plt.figure(figsize = (3.5, 4.5))

    left = 0.15
    height = 0.75
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))
    ax2 = ax.twiny()
    # cax = fig.add_axes([left, bottom+height, width, 0.03])

    parameter, values, label, cmap, logged = variant

    colours = cmr.take_cmap_colors(cmap, len(values), cmap_range=(0.15, 0.85))

    for value, colour in zip(values, colours):

        if logged:
            value_ = 10**value
        else:
            value_ = value

        dagn = relagn(**{parameter: value_})

        Lnu_rel = dagn.get_totSED(rel=True)
        Lnu_nr = dagn.get_totSED(rel=False)

        nu = dagn.nu_obs #frequency grid

        ax.plot(np.log10(nu), np.log10(nu*Lnu_nr), lw='2', color=colour, alpha = 0.3)
        ax.plot(np.log10(nu), np.log10(nu*Lnu_rel), lw=1, color=colour, alpha = 1, label = rf'$\rm {label}={value}$')


    xlims = np.array([14., 20.])
    ylims = (41.1, 45.9)

    ax2.axvline(np.log10(912.), lw=1, c='k', alpha=0.1)

    ax.set_ylim(ylims)
    ax.set_xlim(xlims)

    ax.legend(fontsize=8)

    ax.set_xlabel(r'$\rm \log_{10}(\nu/Hz)$')
    ax.set_ylabel(r'$\rm \log_{10}(\nu L_{\nu}/erg\ s^{-1})$')

    ax2.set_xlim(np.log10(c.to('Angstrom/s').value)-xlims)
    ax2.set_xlabel(r'$\rm \log_{10}(\lambda/\AA)$')

    # add colourbar
    # cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    # fig.colorbar(cmapper, cax=cax, orientation='horizontal')
    # cax.xaxis.tick_top()
    # cax.xaxis.set_label_position('top')
    # cax.set_xlabel(r'$\rm log_{10}(1+\delta_{14})$', fontsize=7)
    # cax.tick_params(axis='x', labelsize=6)


    fig.savefig(f'figs/theory_relagn_{parameter}.pdf')
    fig.clf()