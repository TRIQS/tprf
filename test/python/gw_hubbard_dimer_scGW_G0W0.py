################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by Hugo U.R. Strand
# Author: H. U.R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import sys
import itertools
import numpy as np

from triqs.gf import Gf, MeshImFreq, Idx, inverse
from triqs.gf.gf_factories import make_gf_from_fourier
from triqs.operators import n, c, c_dag, Operator, dagger

from gw_hubbard_dimer import GWHubbardDimer
from gw_hubbard_dimer_matrix import GWHubbardDimerMatrix


def get_ed_g(beta, t, U, wmesh):
    
    from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization
    from triqs.operators import c, c_dag, n, Operator, dagger

    mu = U/2
    H = \
        -mu * (n('up',0) + n('do',0) + n('up',1) + n('do',1) ) + \
        U * n('up',0) * n('do',0) + \
        U * n('up',1) * n('do',1) + \
        -t * ( c_dag('up', 0) * c('up', 1) + c_dag('up', 1) * c('up', 0) ) + \
        -t * ( c_dag('do', 0) * c('do', 1) + c_dag('do', 1) * c('do', 0) )

    fundamental_operators = [ c('up', 0), c('do', 0), c('up', 1), c('do', 1) ]
    ed = TriqsExactDiagonalization(H, fundamental_operators, beta)

    G_w = Gf(mesh=wmesh, target_shape=[2,2])
    G_tau = make_gf_from_fourier(G_w)

    for i, j in itertools.product(range(2), repeat=2):
        ed.set_g2_tau(G_tau[i,j], c('up',i), c_dag('up',j) )
        
    G_w = make_gf_from_fourier(G_tau)

    return G_w


class EDHubbardDimer:

    def __init__(self, beta, t, U, wmesh):

        self.g_w = get_ed_g(beta, t, U, wmesh)
        self.g0_w = get_ed_g(beta, t, 0*U, wmesh)
        self.sigma_w = self.g_w.copy()
        self.sigma_w << inverse(self.g0_w) - inverse(self.g_w)


def compare_gw_solutions(verbose=True):

    beta = 40.0
    U = 0.5
    t = 1.0
    mu = 0.0
    nw = 1024 * 2

    opts = dict(beta=beta, U=U, t=t, mu=mu, nw=nw)

    g0w0     = GWHubbardDimer(maxiter=1, spinless=True,  **opts)
    g0w0_sic = GWHubbardDimer(maxiter=1, **opts)
    
    gw     = GWHubbardDimer(maxiter=100, spinless=True,  **opts)
    gw_sic = GWHubbardDimer(maxiter=100, **opts)
    
    ed = EDHubbardDimer(beta, t, U, gw.g_wk.mesh[0])

    diff_g0w0 = np.max(np.abs(
        g0w0.g_wr[:, Idx(0, 0, 0)][0, 0].data - ed.g_w[0, 0].data))

    diff_scgw = np.max(np.abs(
        gw_sic.g_wr[:, Idx(0, 0, 0)][0, 0].data - ed.g_w[0, 0].data))

    print(f'diff_g0w0 = {diff_g0w0:2.2E}')
    print(f'diff_scgw = {diff_scgw:2.2E}')

    assert( diff_g0w0 < 6e-3 )
    assert( diff_scgw < 8e-5 )
    
    if verbose:
        from triqs.plot.mpl_interface import oplot, oploti, oplotr, plt

        plt.figure(figsize=(8, 5))

        subp = [2, 2, 1]
        xlim = [0, 10]

        plt.subplot(*subp); subp[-1] += 1
        oploti(ed.g_w[0, 0], label='ED')
        oploti(g0w0.g_wr[:, Idx(0, 0, 0)][0, 0], '--', label='G0W0 spin-less')
        oploti(g0w0_sic.g_wr[:, Idx(0, 0, 0)][0, 0], '--', label='G0W0 spin-tensor')
        oploti(gw.g0_wr[:, Idx(0, 0, 0)][0, 0], '-', label='G_0')
        plt.xlim(xlim)
        plt.ylim(top=0)
        plt.ylabel(r'$G(r=0)$')

        plt.subplot(*subp); subp[-1] += 1
        oploti(ed.g_w[0, 0], label='ED')
        oploti(gw.g_wr[:, Idx(0, 0, 0)][0, 0], '--', label='scGW spin-less')
        oploti(gw_sic.g_wr[:, Idx(0, 0, 0)][0, 0], '--', label='scGW spin-tensor')
        oploti(gw.g0_wr[:, Idx(0, 0, 0)][0, 0], '-', label='G_0')
        plt.xlim(xlim)
        plt.ylim(top=0)
        plt.ylabel(r'$G(r=0)$')

        plt.subplot(*subp); subp[-1] += 1
        oploti(ed.sigma_w[0, 0], label='ED')
        oploti(g0w0.sigma_wr[:, Idx(0, 0, 0)][0, 0], '--', label='G0W0 spin-less')
        oploti(g0w0_sic.sigma_wr[:, Idx(0, 0, 0)][0, 0], '--', label='G0W0 spin-tensor')
        plt.xlim(xlim)
        plt.ylim(top=0)
        plt.ylabel(r'$\Sigma(r=0)$')

        plt.subplot(*subp); subp[-1] += 1
        oploti(ed.sigma_w[0, 0], label='ED')
        oploti(gw.sigma_wr[:, Idx(0, 0, 0)][0, 0], '--', label='scGW spin-sum')
        oploti(gw_sic.sigma_wr[:, Idx(0, 0, 0)][0, 0], '--', label='scGW spin-tensor')
        plt.xlim(xlim)
        plt.ylim(top=0)
        plt.ylabel(r'$\Sigma(r=0)$')

        plt.tight_layout()
        #plt.savefig('figure_gw_hubbard_dimer_cf_G0W0_and_scGW_SIC_vs_spinsum.pdf')
        plt.show()

        
if __name__ == '__main__':

    compare_gw_solutions(verbose=False)
