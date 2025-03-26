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

"""
TPSC for the Hubbard model on the square lattice

checking that the TPSC result does satisfy the exact sum rule
Sum_k (chi_sp(k) + chi_ch(k)) = 2n - n^2
at default parameters (U=2, beta=2.5, n=1, square lattice with nearest-neighbour hopping only)

https://jp1.journaldephysique.org/articles/jp1/abs/1997/11/jp1v7p1309/jp1v7p1309.html
"""

import numpy as np

from triqs.gf import MeshDLRImFreq
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.tpsc_solver import tpsc_solver


def test_chi_sum_rule():
    
    # define the parameters
    n = 1.0
    U = 2.0
    wmesh = MeshDLRImFreq(beta=2.5, statistic='Fermion', w_max=12.0, eps=1e-14)
    t = 1.0
    n_k = 128
    # lattice geometry
    units = [(1,0,0), (0,1,0)]
    hoppings = {(+1,+0) : [[-t]],
                (-1,+0) : [[-t]],
                (+0,+1) : [[-t]],
                (+0,-1) : [[-t]]}
    Lat = TBLattice(units=units, hoppings=hoppings)
    kmesh = Lat.get_kmesh(n_k=n_k)
    e_k = Lat.fourier(kmesh)

    # initialize and run solver
    S = tpsc_solver(n=n, U=U, wmesh=wmesh, e_k=e_k, verbose=False)
    S.solve(calc_sigma=False, calc_g=False, check_self_consistency=False)

    # get chi sum rule
    chi_sum = S._get_density(S.chisp_wk + S.chich_wk)

    # test against exact result
    np.testing.assert_array_almost_equal(chi_sum, 2*S.n - S.n**2)

    
if __name__ == '__main__':
    test_chi_sum_rule()
    
