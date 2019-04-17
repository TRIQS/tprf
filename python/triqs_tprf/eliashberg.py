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

from pytriqs.gf import Gf
from lattice import eliashberg_product
from lattice import eliashberg_product_fft
from lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr

# ----------------------------------------------------------------------

def semi_random_initial_delta(g_wk, nr_factor=0.5, seed=None):

    r"""Create a delta based on the GF with random elements

    Returns an anomalous self-energy that can be used as an inital input for the iterative
    solvers. The momentum space is random, while the Matsubara space is only partialy 
    randomized to ensure working tail fits for the Fourier transformations.

    Parameters
    ----------
    g_wk : Gf, 
           Green's function :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the Gf must
           be a MeshProduct with the components (MeshImFreq, MeshBrillouinZone).
    nr_factor : float, optional
                Percentage of :math:`\omega` points which shall not be randomized. This is needed
                to assure a working tail fit for the Fourier transformations. The default is 0.5,
                meaning that 50% of the :math:`\omega` points will not be randomized.
    seed : int, optional
           Set a np.random.seed to enforce predictable results.

    Returns
    -------
    delta : Gf,
            An initial anomalous self-energy :math:`\Delta(i\nu_n, \mathbf{k})` to start
            an iterative solver, given as a Gf with MeshProduct with the components
            (MeshImFreq, MeshBrillouinZone).
    """

    np.random.seed(seed)

    delta = g_wk.copy()
    shape = delta.data.shape
    delta.data[:] = delta.data.real # Pure real delta is sufficient w/o magnetic field
    random_data = np.random.random(shape[1:])
    freq_data = np.mean(np.abs(delta.data), axis=tuple(range(len(shape))[1:]))
    not_randomized = int(nr_factor*shape[0] / 2.)
    start, stop = not_randomized, shape[0]-not_randomized
    freq_data[start:stop] *= np.random.random(stop-start)

    delta.data[:] = np.tensordot(freq_data, random_data, axes=0)

    return delta

def preprocess_gamma_for_fft(Gamma_pp_wk, Gamma_pp_const_k):
    r""" Prepare Gamma to be used with the FFT implementation

    Parameters
    ----------
    Gamma_pp_wk : Gf,
               Pairing vertex :math:`\Gamma(i\omega_n, \mathbf{k})`. The mesh attribute of
               the Gf must be a MeshProduct with the components (MeshImFreq, MeshBrillouinZone).
    Gamma_pp_const_k : float or np.ndarray or Gf
                       Part of the pairing vertex that is constant in Matsubara frequency space
                       :math:`\Gamma(\mathbf{k})`. If given as a Gf its mesh attribute needs to
                       be a MeshBrillouinZone.

    Returns
    -------
    Gamma_pp_dyn_tr : Gf,
                      The dynamic part of Gamma, which converges to zero for
                      :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space.
                      Its mesh attribute is MeshProduct with the components
                      (MeshImTime, MeshCyclicLattice).
    Gamma_pp_const_r : Gf,
                       The constant part of Gamma with mesh attribute MeshCyclicLattice.
    """

    # -- Determine the dynamic and constant part via a tail fit
    # -- (This is done even if the constant term is given to get the specific Gf types)
    Gamma_pp_dyn_wk_fit, Gamma_pp_const_k_fit = split_into_dynamic_wk_and_constant_k(Gamma_pp_wk)

    # -- Use a constant term if explicitly given
    const_type = type(Gamma_pp_const_k)
    if (const_type == float) or (const_type == np.ndarray):
        Gamma_pp_const_k_fit.data[:] = Gamma_pp_const_k
        Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k
    elif (const_type == Gf):
        Gamma_pp_const_k_fit[:] = Gamma_pp_const_k.data
        Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k.data

    # -- FFT dynamic and constant term to (tau, real) or (real)
    Gamma_pp_dyn_tr, Gamma_pp_const_r = dynamic_and_constant_to_tr(Gamma_pp_dyn_wk_fit, 
                                                                    Gamma_pp_const_k_fit)

    return Gamma_pp_dyn_tr, Gamma_pp_const_r

def solve_eliashberg(Gamma_pp_wk, g_wk, initial_delta=None, Gamma_pp_const_k=None, tol=1e-10,
                            product='FFT', solver='PM'):

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
    initial_delta : Gf, optional
                   An initial anomalous self-energy :math:`\Delta(i\nu_n, \mathbf{k})` to start
                   an iterative solver, given as a Gf with MeshProduct with the components
                   (MeshImFreq, MeshBrillouinZone).
                   If not given :func:`semi_random_initial_delta` will be called.
    Gamma_pp_const_k : float or np.ndarray or Gf, optional
                       Part of the pairing vertex that is constant in Matsubara frequency space
                       :math:`\Gamma(\mathbf{k})`. If given as a Gf its mesh attribute needs to
                       be a MeshBrillouinZone.
    tol : float, optional
          Relative accuracy for eigenvalues (stopping criterion).
    product : str, ['FFT', 'SUM'], optional
              Which function of the Eliashberg product shall be used:

                  'FFT' : triqs_tprf.lattice.eliashberg_product_fft,
                          which uses Fourier transformation for optimal computational efficiency.
                          Restrictions : 
                  'SUM' : triqs_tprf.lattice.eliashberg_product, uses the explicit sum.
                          Restrictions : wmesh of Gamma_pp_wk must be atleast twice the size
                                         of the one of g_wk.
    solver : str, ['PM', 'IRAM'], optional
             Which eigenvalue solver shall be used:

                 'PM' : Use the Power Method implemented in :func:`power_method_LR`.

                 'IRAM' : Use the Implicitly Restarted Arnoldi Method implemented in :func:`implicitly_restarted_arnoldi_method`.

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
    :ref:`eliashberg` : Theory of the linearized Eliashberg equation.
    """

    def from_x_to_wk(delta_x):
        delta_wk = g_wk.copy()
        delta_wk.data[:] = delta_x.reshape(delta_wk.data.shape)
        return delta_wk

    def from_wk_to_x(delta_wk):
        delta_x = delta_wk.data.copy().flatten()
        return delta_x

    if product == 'FFT':

        Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(Gamma_pp_wk, Gamma_pp_const_k)

        def matvec(delta_x):
            delta_wk = from_x_to_wk(delta_x)
            delta_out_wk = eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk)
            delta_out_x = from_wk_to_x(delta_out_wk)
            return delta_out_x

    elif product == 'SUM':

        def matvec(delta_x):
            delta_wk = from_x_to_wk(delta_x)
            delta_out_wk = eliashberg_product(Gamma_pp_wk, g_wk, delta_wk)
            delta_out_x = from_wk_to_x(delta_out_wk)
            return delta_out_x

    else:
        raise NotImplementedError('There is no implementation of the eliashberg product'
                                    ' called %s.'%product)
    if not initial_delta:
        initial_delta = semi_random_initial_delta(g_wk)
    initial_delta = from_wk_to_x(initial_delta)

    if solver == 'PM':
        es, evs = power_method_LR(matvec, initial_delta, tol=tol)
        es, evs = [es], [evs]

    elif solver == 'IRAM':
        es, evs = implicitly_restarted_arnoldi_method(matvec, initial_delta, tol=tol)

    else:
        raise NotImplementedError('There is no solver called %s.'%solver)

    eigen_modes = [from_x_to_wk(ele) for ele in evs]

    return es, eigen_modes

def implicitly_restarted_arnoldi_method(matvec, init, tol=1e-10):

    """Find the eigenvalue with the largest real value via the Implicitly Restarted 
    Arnoldi Method

    Parameters
    ----------
    matvec : callable f(v),
             Returns A*v.
    init : np.ndarray,
           The array representation of the anomalous self-energy to start the iterative
           method with. Restriction: len(init.shape) == 1.
    tol : float, optional
          The tolerance at which the iterative scheme is considered to be converged.

    Returns
    -------
    Es : list of float,
           The eigenvalues with the largest positive real part.
    U : list of np.ndarray,
          The corresponding eigenvectors.

    Notes
    -----
    `scipy.sparse.linalg.eigs <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.eigs.html>`_

    `scipy.sparse.linalg.LinearOperator <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html>`_
    """
    N = init.shape[0]
    linop = LinearOperator(matvec=matvec, dtype=np.complex, shape=(N, N))
    Es, U = eigs(linop, which='LR', tol=tol, v0=init)
    Es = Es.real
    
    return list(Es), list(U.T)

def power_method_LR(matvec, init, tol=1e-10, max_it=1e5):

    """Find the eigenvalue with the largest real value via the power method

    Parameters
    ----------
    matvec : callable f(v),
             Returns A*v.
    init : np.ndarray,
           The array representation of the anomalous self-energy to start the iterative
           method with. Restriction: len(init.shape) == 1.
    tol : float, optional
          The tolerance at which the iterative scheme is considered to be converged.
    max_it : float, optional
             The maximum number of iterations that shall be done before a error is raised.

    Returns
    -------
    norm : float,
           The eigenvalue with the largest positive real part.
    v_k : np.ndarray,
          The corresponding eigenvector.
    """

    def iteration(v_k, offset=0.0):
        v_k1 = matvec(v_k) - offset*v_k
        v_k1_norm = np.linalg.norm(v_k1)
        v_k1 = v_k1 / v_k1_norm
        return v_k1_norm+offset, v_k1

    def power_method(init, offset=0.0, tol=tol, max_it=max_it):
        norm, v_k = iteration(init, offset)
        it = 1
        while True:
            norm, new_v_k = iteration(v_k, offset)
            
            # -- Convergence criterion
            add = np.max(np.abs(v_k + new_v_k))
            diff = np.max(np.abs(v_k - new_v_k))

            if (np.allclose(add, 0, atol=tol)) or (np.allclose(diff, 0, atol=tol)):
                break

            v_k = new_v_k
            it += 1
            if it > max_it:
                raise AssertionError('Did not converge.')
        return norm, v_k

    # Find eigenvalue with maximum magnitude
    norm, v_k = power_method(init, tol=tol)

    # Check sign of found eigenvalue
    _, v_k_test = iteration(v_k)

    add = np.sum(np.abs(v_k + v_k_test)) # small if sign of E is negative
    diff = np.sum(np.abs(v_k - v_k_test)) # small if sign of E is positive

    # -- Return eigenvalue with largest real part
    if diff > add: # The eigenvalue with the largest magnitude is negative
        norm, v_k = power_method(init, offset=-norm, tol=tol)
    return norm, v_k

