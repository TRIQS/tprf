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

comparing screened interactions Uch and Usp
with results from the reference implementation
https://github.com/amstremblay/TPSC
"""

import numpy as np

from triqs.mesh import MeshDLRImFreq
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.tpsc_solver import tpsc_solver


def test_tpsc_hubbard_square_lattice():

    # Parameters
    U = 2.0     # Hubbard interaction
    n = 1.0     # Electron density (half-filling)
    t = -1.0    # Nearest neighbour hopping

    wmesh = MeshDLRImFreq(beta=2.5, statistic='Fermion', w_max=12.0, eps=1e-14)
    
    # Square lattice tight binding model
    
    tb = TBLattice(
        units=[(1,0,0), (0,1,0)],
        hoppings={(+1,+0) : [[t]],
                  (-1,+0) : [[t]],
                  (+0,+1) : [[t]],
                  (+0,-1) : [[t]]})

    e_k = tb.fourier(tb.get_kmesh(n_k=16)) # n_k=128 for exact reference
    
    S = tpsc_solver(n, U, wmesh, e_k)
    S.solve()

    # -- Reference results for

    # screened spin and charge interaction vertices

    Usp_ref = 1.5104441009131098
    Uch_ref = 3.597697989091686

    # double occupancy
    
    docc_ref = 0.1888055126141387

    np.testing.assert_array_almost_equal(S.Usp, Usp_ref, decimal=4)
    np.testing.assert_array_almost_equal(S.Uch, Uch_ref, decimal=4)
    np.testing.assert_array_almost_equal(S.docc, docc_ref, decimal=4)


if __name__ == '__main__':
    test_tpsc_hubbard_square_lattice()
