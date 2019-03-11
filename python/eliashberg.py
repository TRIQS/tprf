################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2019 by The Simons Foundation
# Author: H. U.R. Strand
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import eigs

# ----------------------------------------------------------------------

from lattice import eliashberg_product

# ----------------------------------------------------------------------
def solve_eliashberg(Gamma_pp, g_wk, tol=1e-7):

    """ Solve the linearized Eliashberg equation using 
    iterative eigenvalue algorithms from scipy """
    
    def from_x_to_wk(delta_x):
        delta_wk = g_wk.copy()
        delta_wk.data[:] = delta_x.reshape(delta_wk.data.shape)
        return delta_wk

    def from_wk_to_x(delta_wk):
        delta_x = delta_wk.data.copy().flatten()
        return delta_x
    
    def matvec(delta_x):
        delta_wk = from_x_to_wk(delta_x)
        delta_out_wk = eliashberg_product(Gamma_pp, g_wk, delta_wk)
        delta_out_x = from_wk_to_x(delta_out_wk)
        return delta_out_x

    x = from_wk_to_x(g_wk)
    N = x.shape[0]
    linop = LinearOperator(matvec=matvec, dtype=np.complex, shape=(N, N))
    E, U = eigs(linop, which='LR', tol=tol)

    eigen_modes = []
    for idx in xrange(U.shape[-1]):
        delta_wk = from_x_to_wk(U[:, idx])
        eigen_modes.append(delta_wk)

    return E, eigen_modes
    
