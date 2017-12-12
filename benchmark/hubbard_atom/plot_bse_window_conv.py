
""" Plot convergence of the magnetic vertex function in the PH 
channel for the hubbard atom, as a function of fermionic window size,
when computed from numerical inversion of the Bethe Saltpeter equation.

Author: Hugo U.R. Strand (2017), hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np
from pytriqs.gf import *
from analytic_hubbard_atom import analytic_hubbard_atom

# ----------------------------------------------------------------------
if __name__ == '__main__':

    nwf_vec = np.array([5, 10, 20, 40, 80, 160])
    diff_vec = np.zeros_like(nwf_vec, dtype=np.float)

    for idx, nwf in enumerate(nwf_vec):
         d = analytic_hubbard_atom(beta=2.0, U=5.0, nw=1, nwf=nwf)
         diff_vec[idx] = np.max(np.abs( d.gamma_m.data - d.gamma_m_num.data ))
         print 'nwf, diff =', idx, nwf, diff_vec[idx]

    # ------------------------------------------------------------------
    # -- Plot
    
    from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

    x = 1./nwf_vec

    plt.figure(figsize=(3.25, 3))
    
    plt.plot(x, diff_vec, 'o-', alpha=0.75)

    plt.xlabel(r'$1/n_{w}$')
    plt.ylabel(r'$\max|\Gamma_{ana} - \Gamma_{num}|$')

    plt.ylim([0, diff_vec.max()])
    plt.xlim([0, x.max()])

    plt.tight_layout()
    plt.savefig('figure_bse_hubbard_atom_convergence.pdf')

    plt.show()
