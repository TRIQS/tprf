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
Testing the improved bubble functionality necessary for TPSC+
DOI: https://doi.org/10.1103/PhysRevB.108.075144
"""

import numpy as np

from triqs.gf import MeshDLRImFreq
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.tpsc_solver import tpsc_solver

def test_GG0_bubble():

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
    S.solve(calc_sigma=True, calc_g=True, check_self_consistency=False)

    # get improved bubble
    S._imtime_bubble_chi2_wk()