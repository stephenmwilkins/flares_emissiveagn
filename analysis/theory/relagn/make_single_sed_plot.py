
import numpy as np
import matplotlib.pyplot as plt
from relagn import relagn
from unyt import c

# set style
plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')


dagn = relagn(a=0.0, cos_inc=0.5, log_mdot = 0.5, M=1E9) 


Lnu_rel = dagn.get_totSED(rel=True)
Lnu_nr = dagn.get_totSED(rel=False)


fig = plt.figure(figsize = (3.5, 4.5))

left = 0.15
height = 0.75
bottom = 0.15
width = 0.8

ax = fig.add_axes((left, bottom, width, height))
ax2 = ax.twiny()
# cax = fig.add_axes([left, bottom+height, width, 0.03])


nu = dagn.nu_obs #frequency grid

ax.plot(np.log10(nu), np.log10(nu*Lnu_nr), ls='dashed', color='grey')

ax.plot(np.log10(nu), np.log10(nu*Lnu_rel), color='k')


xlims = np.array([14., 20.])
ylims = (45., 48.)

ax.set_ylim(ylims)
ax.set_xlim(xlims)



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


fig.savefig(f'figs/single.pdf')


fig.clf()