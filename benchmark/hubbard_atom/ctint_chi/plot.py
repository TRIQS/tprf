# ----------------------------------------------------------------------

import numpy as np

from pytriqs.gf import *
from h5 import HDFArchive
from pytriqs.utility import mpi

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from itertools import product

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.analytic_hubbard_atom import analytic_hubbard_atom

from pytriqs.plot.mpl_interface import plt, oplot, oploti, oplotr

# ----------------------------------------------------------------------
if __name__ == '__main__':

    with HDFArchive("data_ctint.h5", 'r') as results:
        p = results["p"]

    print p
    p.G_iw['up'].name = r'$G_{ctint, \uparrow}$'
    p.G_iw['dn'].name = r'$G_{ctint, \downarrow}$'

    p.chi_m = p.G2_iw[('up', 'up')] - p.G2_iw[('up', 'dn')]
    
    parm = ParameterCollection(
        U = p.U,
        beta = p.beta,
        nw = 1,
        nwf = p.n_iw / 2,
        nwf_gf = p.n_iw,
        )
    
    a = analytic_hubbard_atom(**parm.dict())
    a.G_iw.name = r'$G_{analytic}$'

    plt.figure(figsize=(3.25*2, 3*2))

    subp = [2, 2, 1]

    plt.subplot(*subp); subp[-1] += 1
    
    oploti(p.G_iw)
    oploti(a.G_iw)
    
    plt.subplot(*subp); subp[-1] += 1
    diff = a.G_iw.copy()
    diff << p.G_iw['up'] - a.G_iw
    diff.name = r'$G_{ctint} - G_{analytic}$'
    oplot(diff)

    plt.subplot(*subp); subp[-1] += 1
    vmax = np.max(np.abs(p.chi_m.data.real))
    opt = dict(vmax=vmax, vmin=-vmax, cmap='PuOr')
    data = np.squeeze(p.chi_m.data.real)
    plt.pcolormesh(data, **opt)
    
    plt.tight_layout()
    plt.show()
