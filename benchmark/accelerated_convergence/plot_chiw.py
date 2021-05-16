from common import *
from triqs.plot.mpl_interface import oplot, oplotr, plt

from triqs.gf import Fourier, inverse,GfImTime,GfImFreq,SemiCircular,iOmega_n,GfReFreq


with HDFArchive('data_chi.h5', 'r') as a: p = a['p']

#%%
beta = p.chi_tau.mesh.beta

plt.figure(figsize=(3.25*2.2, 2))
subp = [1, 2, 1]
plt.subplot(*subp)
oplot(p.chi_tau.real)

chi_w = GfImFreq(target_shape=p.chi_tau.target_shape,beta = beta,statistic='Boson',n_points=10,name='chi_w')
chi_w << Fourier(p.chi_tau)


subp = [1, 2, 2]
plt.subplot(*subp)

oplot(chi_w.real)

plt.tight_layout()
plt.show()
