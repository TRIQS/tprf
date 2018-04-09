import numpy as np

from pytriqs.gf import Gf, MeshImFreq, MeshImTime
from pytriqs.gf import GfImFreq, GfImTime
from pytriqs.gf import iOmega_n, inverse, InverseFourier

beta = 1.234
eps = 3.55
nw = 100
nt = 1000

g_iw = GfImFreq(name=r'$g(i\omega_n)$', beta=beta,
                 statistic='Fermion', n_points=nw,
                 indices=[1])

g_iw << inverse(iOmega_n - eps)

g_tau = GfImTime(name=r'$g(\tau)$', beta=beta,
                 statistic='Fermion', n_points=nt,
                 indices=[1])

g_tau << InverseFourier(g_iw)

# -- Test fourier

from triqs_tprf.fourier import g_iw_from_tau
g_iw_ref = g_iw_from_tau(g_tau, nw)

np.testing.assert_array_almost_equal(g_iw_ref.data, g_iw.data)

if False:
    from pytriqs.plot.mpl_interface import oplot, plt

    subp = [1, 3, 1]
    plt.subplot(*subp); subp[-1] += 1
    oplot(g_iw)

    plt.subplot(*subp); subp[-1] += 1
    oplot(g_tau)

    plt.subplot(*subp); subp[-1] += 1
    oplot(g_iw_ref)
    
    plt.show()
