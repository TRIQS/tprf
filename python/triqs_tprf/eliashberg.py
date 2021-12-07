# -*- coding: utf-8 -*-

################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2019, The Simons Foundation and S. Käser
# Authors: H. U.R. Strand, S. Käser
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

from .lattice import eliashberg_product
from .lattice import eliashberg_product_fft
from .lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr

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
    for idx in range(U.shape[-1]):
        delta_wk = from_x_to_wk(U[:, idx])
        eigen_modes.append(delta_wk)

    return E, eigen_modes
    
def solve_eliashberg_fft(Gamma_pp_wk, g_wk, Gamma_pp_const_k=None, tol=1e-10):

    r""" Solve the linearized Eliashberg equation
    
    Returns the biggest eigenvalues and corresponding eigenvectors of the linearized Eliashberg
    equation. The Eliashberg equation implementation is using fourier transformations for 
    computational efficiency. The eigenvalues are found using an iterative algorithm from scipy.
    
    Parameters
    ----------
    Gamma_pp_wk : Gf,
               Pairing vertex :math:`\Gamma(i\omega_n, \mathbf{k})`. The mesh attribute of
               the Gf must be a MeshProduct with the components (MeshImFreq, MeshBrillouinZone).
    g_wk : Gf, 
           Green's function :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the Gf must
           be a MeshProduct with the components (MeshImFreq, MeshBrillouinZone).
    Gamma_pp_const_k : float or np.ndarray or Gf, optional
                       Part of the pairing vertex that is constant in Matsubara frequency space
                       :math:`\Gamma(\mathbf{k})`. If given as a Gf its mesh attribute needs to
                       be a MeshBrillouinZone.
    tol : float, optional
          Relative accuracy for eigenvalues (stopping criterion).

    Returns
    -------
    Es : list of float,
         Biggest eigenvalues of the linearized Eliashberg equation :math:`\lambda`.
    eigen_modes : list of Gf,
                  Corresponding eigenvectors (anomalous self-energies) 
                  :math:`\Delta(i\nu_n, \mathbf{k})` as Gf with MeshProduct with the components
                  (MeshImFreq, MeshBrillouinZone).

    See Also
    --------
    :ref:`eliashberg`

    """

    # -- Determine the dynamic and constant part via a tail fit
    # -- (This is done even if the constant term is given to get the specific GF types)
    Gamma_pp_dyn_wk_fit, Gamma_pp_const_k_fit = split_into_dynamic_wk_and_constant_k(Gamma_pp_wk)

    # -- Use a constant term if explicitly given
    if Gamma_pp_const_k:
        try:
            Gamma_pp_const_k_fit.data[:] = Gamma_pp_const_k
            Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k
        except TypeError:
            Gamma_pp_const_k_fit[:] = Gamma_pp_const_k.data
            Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k.data

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
    Es, U = eigs(linop, which='LR', tol=tol, v0=v0)

    eigen_modes = []
    for idx in range(U.shape[-1]):
        delta_wk = from_x_to_wk(U[:, idx])
        eigen_modes.append(delta_wk)

    return list(Es), eigen_modes
    
