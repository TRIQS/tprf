################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2025 by Xaver Landerl
# Author: Xaver Landerl
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

r"""
TPSC for the Hubbard model on the square lattice

checking that the TPSC result does satisfy the exact sum rules

Sum_k (chi_sp(k) + chi_ch(k)) = 2n - n^2

and

Sum_k (\Sigma^(2)(k) G^(0)(k)) = U * doubleocc

at default parameters (U=2, beta=2.5, n=1, square lattice with nearest-neighbour hopping only)

https://jp1.journaldephysique.org/articles/jp1/abs/1997/11/jp1v7p1309/jp1v7p1309.html
"""

import numpy as np

from triqs.mesh import MeshDLRImFreq
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.tpsc_solver import tpsc_solver


def test_tpsc_sum_rules():

    n = 1.0
    U = 2.0
    t = 1.0

    wmesh = MeshDLRImFreq(beta=2.5, statistic='Fermion', w_max=12.0, eps=1e-14)

    tb = TBLattice(
        units=[(1,0,0), (0,1,0)],
        hoppings={(+1,+0) : [[-t]],
                  (-1,+0) : [[-t]],
                  (+0,+1) : [[-t]],
                  (+0,-1) : [[-t]]})

    e_k = tb.fourier(tb.get_kmesh(n_k=8))
    
    S = tpsc_solver(n, U, wmesh, e_k)
    S.solve(calc_sigma=True)

    tr_chi_sum = S._get_density(S.chisp_wk + S.chich_wk)
    np.testing.assert_array_almost_equal(tr_chi_sum, 2*S.n - S.n**2)
    
    trace_SigmaG0 = S._get_density(S.sigma_wk * S.g0_wk)
    np.testing.assert_array_almost_equal(trace_SigmaG0, S.U*S.docc)


if __name__ == '__main__':
    test_tpsc_sum_rules()
    
