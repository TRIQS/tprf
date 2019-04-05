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
from lattice import eliashberg_product_fft
from lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr

# ----------------------------------------------------------------------
def solve_eliashberg(Gamma_pp, g_wk, tol=1e-10):

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

    np.random.seed(1337)
    v0 = np.random.random(N)
    E, U = eigs(linop, which='LR', tol=tol, v0=v0)

    eigen_modes = []
    for idx in xrange(U.shape[-1]):
        delta_wk = from_x_to_wk(U[:, idx])
        eigen_modes.append(delta_wk)

    return E, eigen_modes
    
def solve_eliashberg_fft(Gamma_pp, g_wk, Gamma_pp_const_k=None, tol=1e-10):

    """ Solve the linearized Eliashberg equation using 
    iterative eigenvalue algorithms from scipy via the FFT product
    
    Parmeters:

        Gamma_pp : gf object, pairing vertex
        g_wk : gf object, Green's function
        Gamma_pp_const_k : int, float, np.ndarray, part of the pairing vertex that is constant
                            in Matsubara frequency space
        tol : float, relative accuracy for eigenvalues (stopping criterion)
    """

    # -- Determine the dynamic and constant part via a tail fit
    # -- (This is done even if the constant term is given to get the specific GF types)
    Gamma_pp_dyn_wk_fit, Gamma_pp_const_k_fit = split_into_dynamic_wk_and_constant_k(Gamma_pp)

    # -- Use a the constant term if explicitly given
    if Gamma_pp_const_k:
        Gamma_pp_const_k_fit.data[:] = Gamma_pp_const_k
        Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp.data - Gamma_pp_const_k

    # -- FFT dynamic and constant term to (tau, real) or (real)
    Gamma_pp_dyn_tr, Gamma_pp_const_r = dynamic_and_constant_to_tr(Gamma_pp_dyn_wk_fit, 
                                                                    Gamma_pp_const_k_fit)
    
    def from_x_to_wk(delta_x):
        delta_wk = g_wk.copy()
        delta_wk.data[:] = delta_x.reshape(delta_wk.data.shape)
        return delta_wk

    def from_wk_to_x(delta_wk):
        delta_x = delta_wk.data.copy().flatten()
        return delta_x
    
    def matvec(delta_x):
        delta_wk = from_x_to_wk(delta_x)
        delta_out_wk = eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk)
        delta_out_x = from_wk_to_x(delta_out_wk)
        return delta_out_x

    x = from_wk_to_x(g_wk)
    N = x.shape[0]
    linop = LinearOperator(matvec=matvec, dtype=np.complex, shape=(N, N))

    np.random.seed(1337)
    v0 = np.random.random(N)
    E, U = eigs(linop, which='LR', tol=tol, v0=v0)

    eigen_modes = []
    for idx in xrange(U.shape[-1]):
        delta_wk = from_x_to_wk(U[:, idx])
        eigen_modes.append(delta_wk)

    return E, eigen_modes
    
