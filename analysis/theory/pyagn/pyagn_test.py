from pyagn.sed import SED
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from unyt import c, Angstrom, Hz, h, keV


# initialize class
M = 1e8
mdot = 0.5
a = 0
sed_test = SED(M=M, mdot=mdot, astar=0)

sed = sed_test.disk_sed


# lam = np.arange(1, 10000, 1) * Angstrom
# energy = h*c/lam


log10energy = np.arange(-5., 3., 0.1)
energy = 10**log10energy * keV

luminosity = np.array([sed_test.disk_spectral_luminosity(e) for e in energy.to('erg').value])


plt.plot(log10energy, np.log10(luminosity))
plt.xlim([-5, 3])
plt.show()