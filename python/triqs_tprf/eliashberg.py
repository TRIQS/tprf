# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019, The Simons Foundation and S. Käser
# Author: S. Käser, H. U.R. Strand
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

import functools
import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import eigs

# ----------------------------------------------------------------------

from triqs.gf import Gf
from triqs.gf.meshes import MeshDLRImFreq
from .lattice import eliashberg_product
from .lattice import eliashberg_product_fft, eliashberg_product_fft_constant
from .lattice import split_into_dynamic_wk_and_constant_k, dynamic_and_constant_to_tr
from .lattice import construct_phi_wk

# ----------------------------------------------------------------------


def solve_eliashberg(
    Gamma_pp_wk,
    g_wk,
    initial_delta=None,
    Gamma_pp_const_k=None,
    tol=1e-10,
    product="FFT",
    solver="IRAM",
    symmetrize_fct=lambda x: x,
    k=6,
    fmp_index = 0,
):
    r""" Solve the linearized Eliashberg equation
    
    Returns the biggest eigenvalues and corresponding eigenvectors of the linearized Eliashberg
    equation, for a particle-particle vertex in the random phase approximation,
    as described here :ref:`eliashberg_rpa`. The Eliashberg equation implementation is
    using fourier transformations for computational efficiency. The eigenvalues are found
    using an iterative algorithm from scipy.
    
    Parameters
    ----------
    Gamma_pp_wk : Gf,
               Pairing vertex :math:`\Gamma(i\omega_n, \mathbf{k})`. The mesh attribute of
               the Gf must be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    g_wk : Gf, 
           Green's function :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the Gf must
           be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    initial_delta : Gf, optional
                   An initial anomalous self-energy :math:`\Delta(i\nu_n, \mathbf{k})` to start
                   an iterative solver, given as a Gf with MeshProduct with the components
                   (MeshImFreq, MeshBrZone).
                   If not given :func:`semi_random_initial_delta` will be called.
    Gamma_pp_const_k : float or np.ndarray or Gf, optional
                       Part of the pairing vertex that is constant in Matsubara frequency space
                       :math:`\Gamma(\mathbf{k})`. If given as a Gf its mesh attribute needs to
                       be a MeshBrZone. If not given, the constant part will be fitted.
    tol : float, optional
          Relative accuracy for eigenvalues (stopping criterion).
    product : str, ['FFT', 'SUM'], optional
              Which function of the Eliashberg product shall be used:

                  'FFT' : triqs_tprf.lattice.eliashberg_product_fft,
                          which uses Fourier transformation for optimal computational efficiency.

                  'SUM' : triqs_tprf.lattice.eliashberg_product, uses the explicit sum.
                          Restrictions : wmesh of Gamma_pp_wk must be atleast twice the size of the one of g_wk.

    solver : str, ['IRAM', 'PM'], optional
             Which eigenvalue solver shall be used:

                 'IRAM' : Use the Implicitly Restarted Arnoldi Method implemented in :func:`implicitly_restarted_arnoldi_method`.

                 'PM' : Use the Power Method implemented in :func:`power_method_LR`.

    symmetrize_fct : function, optional
                     A function that takes one parameter: A Green's function 
                     :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the
                     Gf must be a MeshProduct with the components 
                     (MeshImFreq, MeshBrZone).
                     This function is applied after every iteration of the
                     eigenvalue solver and can be used to enforce a specific
                     symmetry. If no symmetries are enforced, caution is need, because
                     unphysical symmetries can occur.

    k : int, optional
        The number of leading superconducting gaps that shall be calculated. Does
        only have an effect, if 'IRAM' is used as a solver.

    fmp_index : int, optional
                The momentum mesh data index corresponding to the momentum used in a 
                finite momentum pairing calculation. The default value corresponds to no
                finite momentum pairing.

    Returns
    -------
    Es : list of float,
         Biggest eigenvalues of the linearized Eliashberg equation :math:`\lambda`.
    eigen_modes : list of Gf,
                  Corresponding eigenvectors (anomalous self-energies) 
                  :math:`\Delta(i\nu_n, \mathbf{k})` as Gf with MeshProduct with the components
                  (MeshImFreq, MeshBrZone).

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

    hasDLRMesh = type(Gamma_pp_wk.mesh.components[0]) == MeshDLRImFreq

    if product == "FFT":

        Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(
            Gamma_pp_wk, Gamma_pp_const_k
        )

        if np.allclose(
            Gamma_pp_dyn_tr.data, 0
        ):  # -- If dynamic part is zero reduced calculation
            eli_prod = functools.partial(
                eliashberg_product_fft_constant, Gamma_pp_const_r, g_wk, fmpindex=fmp_index
            )

        else:
            eli_prod = functools.partial(
                eliashberg_product_fft, Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, fmpindex=fmp_index
            )

    elif product == "SUM":
        if(hasDLRMesh): 
            raise NotImplementedError(
                "There is no implementation of the eliashberg product "
                "called %s when using DLR. Please use the FFT product instead."%product
            )

        eli_prod = functools.partial(eliashberg_product, Gamma_pp_wk, g_wk, fmpindex=fmp_index)

    else:
        raise NotImplementedError(
            "There is no implementation of the eliashberg product"
            " called %s." % product
        )

    def matvec(delta_x):
        delta_wk = from_x_to_wk(delta_x)
        delta_out_wk = eli_prod(delta_wk)
        delta_out_wk = symmetrize_fct(delta_out_wk)
        delta_out_x = from_wk_to_x(delta_out_wk)
        return delta_out_x

    if not initial_delta:
        initial_delta = semi_random_initial_delta(g_wk)
    initial_delta = from_wk_to_x(initial_delta)

    if solver == "PM":
        es, evs = power_method_LR(matvec, initial_delta, tol=tol)
        es, evs = [es], [evs]

    elif solver == "IRAM":
        es, evs = implicitly_restarted_arnoldi_method(
            matvec, initial_delta, k=k, tol=tol
        )

    else:
        raise NotImplementedError("There is no solver called %s." % solver)

    eigen_modes = [from_x_to_wk(ele) for ele in evs]

    return es, eigen_modes


def preprocess_gamma_for_fft(Gamma_pp_wk, Gamma_pp_const_k=None):
    r""" Prepare Gamma to be used with the FFT implementation

    Parameters
    ----------
    Gamma_pp_wk : Gf,
               Pairing vertex :math:`\Gamma(i\omega_n, \mathbf{k})`. The mesh attribute of
               the Gf must be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    Gamma_pp_const_k : float or np.ndarray or Gf
                       Part of the pairing vertex that is constant in Matsubara frequency space
                       :math:`\Gamma(\mathbf{k})`. If given as a Gf its mesh attribute needs to
                       be a MeshBrZone. If not given, the constant part will be fitted.

    Returns
    -------
    Gamma_pp_dyn_tr : Gf,
                      The dynamic part of Gamma, which converges to zero for
                      :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space.
                      Its mesh attribute is MeshProduct with the components
                      (MeshImTime, MeshCycLat).
    Gamma_pp_const_r : Gf,
                       The constant part of Gamma with mesh attribute MeshCycLat.
    """

    # -- If the function has a DLR mesh we cannot use split_into_dynamic_wk_and_constant_k yet
    # TODO: implement split_into_dynamic_wk_and_constant_k for DLR-mesh Gfs
    hasDLRMesh = type(Gamma_pp_wk.mesh.components[0]) == MeshDLRImFreq
    if(hasDLRMesh):
        if(Gamma_pp_const_k is None):
            Gamma_pp_const_k = Gf(mesh=Gamma_pp_wk.mesh.components[1], target_shape=Gamma_pp_wk.target_shape)
            Gamma_pp_const_k.data[:] = 0.0

        Gamma_pp_dyn_wk = Gamma_pp_wk.copy()
        Gamma_pp_dyn_wk.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k.data
        Gamma_pp_dyn_tr, Gamma_pp_const_r = dynamic_and_constant_to_tr(
            Gamma_pp_dyn_wk, Gamma_pp_const_k
        )

        return Gamma_pp_dyn_tr, Gamma_pp_const_r

    # -- Determine the dynamic and constant part via a tail fit
    # -- (This is done even if the constant term is given to get the specific Gf types)
    Gamma_pp_dyn_wk_fit, Gamma_pp_const_k_fit = split_into_dynamic_wk_and_constant_k(
        Gamma_pp_wk
    )

    # -- Use a constant term if explicitly given
    const_type = type(Gamma_pp_const_k)
    if (const_type == float) or (const_type == np.ndarray):
        Gamma_pp_const_k_fit.data[:] = Gamma_pp_const_k
        Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k
    elif const_type == Gf:
        Gamma_pp_const_k_fit.data[:] = Gamma_pp_const_k.data
        Gamma_pp_dyn_wk_fit.data[:] = Gamma_pp_wk.data - Gamma_pp_const_k.data
    # -- FFT dynamic and constant term to (tau, real) or (real)
    Gamma_pp_dyn_tr, Gamma_pp_const_r = dynamic_and_constant_to_tr(
        Gamma_pp_dyn_wk_fit, Gamma_pp_const_k_fit
    )

    return Gamma_pp_dyn_tr, Gamma_pp_const_r


def semi_random_initial_delta(g_wk, nr_factor=0.5, seed=None):
    r"""Create a delta based on the GF with random elements

    Returns an anomalous self-energy that can be used as an inital input for the iterative
    solvers. The momentum space is random, while the Matsubara space is only partialy 
    randomized to ensure working tail fits for the Fourier transformations.

    Parameters
    ----------
    g_wk : Gf, 
           Green's function :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the Gf must
           be a MeshProduct with the components (MeshImFreq, MeshBrZone).
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
            (MeshImFreq, MeshBrZone).
    """

    np.random.seed(seed)

    delta = g_wk.copy()
    shape = delta.data.shape
    delta.data[:] = delta.data.real  # Pure real delta is sufficient w/o magnetic field
    random_data = np.random.random(shape[1:])
    freq_data = np.mean(np.abs(delta.data), axis=tuple(range(len(shape))[1:]))
    not_randomized = int(nr_factor * shape[0] / 2.0)
    start, stop = not_randomized, shape[0] - not_randomized
    freq_data[start:stop] *= np.random.random(stop - start)

    delta.data[:] = np.tensordot(freq_data, random_data, axes=0)

    return delta


def implicitly_restarted_arnoldi_method(matvec, init, tol=1e-10, k=6):
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
    k : int, optional
        The number of eigenvalues and eigenvectors desired.

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
    linop = LinearOperator(matvec=matvec, dtype=complex, shape=(N, N))
    Es, U = eigs(linop, k=k, which="LR", tol=tol, v0=init)
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
        v_k1 = matvec(v_k) - offset * v_k
        v_k1_norm = np.linalg.norm(v_k1)
        v_k1 = v_k1 / v_k1_norm
        return v_k1_norm + offset, v_k1

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
                raise AssertionError("Did not converge.")
        return norm, v_k

    # Find eigenvalue with maximum magnitude
    norm, v_k = power_method(init, tol=tol)

    # Check sign of found eigenvalue
    _, v_k_test = iteration(v_k)

    add = np.sum(np.abs(v_k + v_k_test))  # small if sign of E is negative
    diff = np.sum(np.abs(v_k - v_k_test))  # small if sign of E is positive

    # -- Return eigenvalue with largest real part
    if diff > add:  # The eigenvalue with the largest magnitude is negative
        norm, v_k = power_method(init, offset=-norm, tol=tol)
    return norm, v_k


def allclose_by_scalar_multiplication(delta_1, delta_2, atol=1e-10):
    """Test if two eigenvectors are equal if multiplied by a scalar

    Eigenvectors are not unique and can be multiplied by any complex scalar.
    Therfore two eigenvalue solver could output different eigenvectors for 
    the same non-degenerate eigenvalue.
    This function checks if two eigenvectors are only different, because of a multiplication
    by a scalar.

    Parameters
    ----------
    delta_1 : Gf
    delta_2 : Gf
    tol : float, optional
          The tolerance at which the eigenvector are considered to be equal up to a scalar 

    Returns
    -------
    have_common_scalar_factor : bool,
                                True if the two eigenvectors are equal up to a scalar.
                                False otherwise.
    """
    delta_1_arr = delta_1.data.flatten()
    delta_2_arr = delta_2.data.flatten()

    # Remove numerical zeroes
    delta_1_arr = delta_1_arr[np.abs(delta_1_arr) > 1e-7]
    delta_2_arr = delta_2_arr[np.abs(delta_2_arr) > 1e-7]

    try:
        division_of_deltas = np.divide(delta_1_arr, delta_2_arr)
    except ValueError:  # Arrays do not contain the same # of zeroes and are therefore not equal
        return False

    # Check if elements share common scalar factor
    have_common_scalar_factor = np.allclose(
        division_of_deltas, division_of_deltas[0], atol=atol
    )

    return have_common_scalar_factor


def construct_gamma_singlet_rpa(U_d, U_m, phi_d_wk, phi_m_wk):
    r"""Construct the irreducible singlet vertex in the RPA limit

    The irreducible singlet vertex in the random phase approximation limit for a
    symmetrized calculations of the Eliashberg equation is given by

    .. math::
        \Gamma^{\text{s}}_{a\overline{b}c\overline{d}}(Q=0, K, K') \equiv
        \frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{d}}
        +
        \frac{3}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{m}}
        +
        \text{Complex Conjugate}
        \left[
        3 
        \Phi^{\text{m}}_{c\overline{b}a\overline{d}}(K-K')
        -
        \Phi^{\text{d}}_{c\overline{b}a\overline{d}}(K-K')
        \right]
        \,.

    For more details see :ref:`eliashberg`.

    Parameters
    ----------
    U_d : np.ndarray,
          The local static interaction in the density channel.
    U_m : np.ndarray,
          The local static interaction in the magnetic channel.
    phi_d_wk : Gf,
               The reducible ladder vertex in the density channel
               :math:`\Phi^{\mathrm{d}}(i\omega_n, \mathbf{q})`. The mesh attribute of the Gf
               must be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    phi_m_wk : Gf,
               The reducible ladder vertex in the magnetic channel
               :math:`\Phi^{\mathrm{m}}(i\omega_n, \mathbf{q})`. The mesh attribute of the Gf
               must be a MeshProduct with the components (MeshImFreq, MeshBrZone).

    Returns
    -------
    gamma_singlet : Gf,
                    The irreducible singlet vertex in the RPA limit for a symmetrized
                    calculation of the Eliashberg equation
                    :math:`\Gamma^{\mathrm{s}}(i\omega_n,\mathbf{q})`.
    """
    gamma_singlet = 0.0 * phi_d_wk.copy()

    gamma_singlet.data[:] = 3 * np.conjugate(phi_m_wk.data) + np.conjugate(phi_d_wk.data)
    gamma_singlet.data[:] = gamma_singlet.data.transpose([0, 1, 4, 3, 2, 5])

    gamma_singlet.data[:] += 0.5 * U_d + 1.5 * U_m
    return gamma_singlet


def construct_gamma_triplet_rpa(U_d, U_m, phi_d_wk, phi_m_wk):
    r"""Construct the irreducible triplet vertex in the RPA limit

    The irreducible triplet vertex in the random phase approximation limit for a
    symmetrized calculations of the Eliashberg equation is given by

    .. math::
        \Gamma^{\text{t}}_{a\overline{b}c\overline{d}}(Q=0, K, K') \equiv
        -
        \frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{d}}
        +
        \frac{1}{2}U_{a\overline{b}c\overline{d}}^{\mathrm{m}} 
        +
        \text{Complex Conjuate}
        \left[
        -
        \Phi^{\text{m}}_{c\overline{b}a\overline{d}}(K-K')
        -
        \Phi^{\text{d}}_{c\overline{b}a\overline{d}}(K-K')
        \right]
        \,.

    For more details see :ref:`eliashberg`.

    Parameters
    ----------
    U_d : np.ndarray,
          The local static interaction in the density channel.
    U_m : np.ndarray,
          The local static interaction in the magnetic channel.
    phi_d_wk : Gf,
               The reducible ladder vertex in the density channel
               :math:`\Phi^{\mathrm{d}}(i\omega_n, \mathbf{q})`. The mesh attribute of the Gf
               must be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    phi_m_wk : Gf,
               The reducible ladder vertex in the magnetic channel
               :math:`\Phi^{\mathrm{m}}(i\omega_n, \mathbf{q})`. The mesh attribute of the Gf
               must be a MeshProduct with the components (MeshImFreq, MeshBrZone).

    Returns
    -------
    gamma_triplet : Gf,
                    The irreducible triplet vertex in the RPA limit for a symmetrized
                    calculation of the Eliashberg equation
                    :math:`\Gamma^{\mathrm{t}}(i\omega_n,\mathbf{q})`.
    """
    gamma_triplet = 0.0 * phi_d_wk.copy()

    gamma_triplet.data[:] = -np.conjugate(phi_m_wk.data) - np.conjugate(phi_d_wk.data)
    gamma_triplet.data[:] = gamma_triplet.data.transpose([0, 1, 4, 3, 2, 5])

    gamma_triplet.data[:] += -0.5 * U_d + 0.5 * U_m
    return gamma_triplet



def solve_eliashberg_nonlinear(
    Gamma_pp_wk,
    g_wk,
    initial_delta=None,
    Gamma_pp_const_k=None,
    tol=1e-10,
    max_it=1e5,
    mixing_frac = 0.1,
    steffensens_method = False,
    product="FFT",
    symmetrize_fct=lambda x: x,
    fmp_index = 0,
    verbose = False,
):
    r""" Solve the non-linearized single-band Eliashberg equation
   
    Returns the anomalous self-energy, obtained by solving the non-linearized single-band
    Eliashberg equation. The Eliashberg equation implementation is using fourier transformations
    for computational efficiency. The anomalous self-energy is found using an iterative
    root-finding algorithm.

    Parameters
    ----------
    Gamma_pp_wk : Gf,
               Pairing vertex :math:`\Gamma(i\omega_n, \mathbf{k})`. The mesh attribute of
               the Gf must be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    g_wk : Gf, 
           Green's function :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the Gf must
           be a MeshProduct with the components (MeshImFreq, MeshBrZone).
    initial_delta : Gf, optional
                   An initial anomalous self-energy :math:`\Delta(i\nu_n, \mathbf{k})` to start
                   an iterative solver, given as a Gf with MeshProduct with the components
                   (MeshImFreq, MeshBrZone).
                   If not given :func:`semi_random_initial_delta` will be called.
    Gamma_pp_const_k : float or np.ndarray or Gf, optional
                       Part of the pairing vertex that is constant in Matsubara frequency space
                       :math:`\Gamma(\mathbf{k})`. If given as a Gf its mesh attribute needs to
                       be a MeshBrZone. If not given, the constant part will be fitted.
    tol : float, optional
          Relative accuracy for anomalous self-energy (stopping criterion).
    max_it : int, optional
          Maximum number of iterations performed.
    mixing_frac : float, optional
          Fraction of the delta of the previous iteration is mixed in with the new delta to
          avoid trapping of the algorithm.
    steffensens_method : bool, optional
          Whether to use Steffensen's method as the iterative root-finding algorithm.
    product : str, ['FFT', 'SUM'], optional
              Which function of the Eliashberg product shall be used:

                  'FFT' : triqs_tprf.lattice.eliashberg_product_fft,
                          which uses Fourier transformation for optimal computational efficiency.

                  'SUM' : triqs_tprf.lattice.eliashberg_product, uses the explicit sum.
                          Restrictions : wmesh of Gamma_pp_wk must be atleast twice the size of the one of g_wk.
    symmetrize_fct : function, optional
                     A function that takes one parameter: A Green's function 
                     :math:`G(i\nu_n, \mathbf{k})`. The mesh attribute of the
                     Gf must be a MeshProduct with the components 
                     (MeshImFreq, MeshBrZone).
                     This function is applied after every iteration of the
                     solver and can be used to enforce a specific symmetry.
                     If no symmetries are enforced, caution is need, because
                     unphysical symmetries can occur.
    fmp_index : int, optional
                The momentum mesh data index corresponding to the momentum used in a 
                finite momentum pairing calculation. The default value corresponds to no
                finite momentum pairing.

    Returns
    -------
    delta_wk : Gf,
               The anomalous self-energy :math:`\Delta(i\nu_n, \mathbf{k})` 
               as Gf with MeshProduct with the components (MeshImFreq, MeshBrZone).
    """
    
    def from_x_to_wk(delta_x):
        delta_wk = g_wk.copy()
        delta_wk.data[:] = delta_x.reshape(delta_wk.data.shape)
        return delta_wk

    def from_wk_to_x(delta_wk):
        delta_x = delta_wk.data.copy().flatten()
        return delta_x

    hasDLRMesh = type(Gamma_pp_wk.mesh.components[0]) == MeshDLRImFreq

    if product == "FFT":

        Gamma_pp_dyn_tr, Gamma_pp_const_r = preprocess_gamma_for_fft(
            Gamma_pp_wk, Gamma_pp_const_k
        )

        if np.allclose(
            Gamma_pp_dyn_tr.data, 0
        ):  # -- If dynamic part is zero reduced calculation
            eli_prod = lambda delta_wk: eliashberg_product_fft_constant(Gamma_pp_const_r, g_wk, delta_wk, linearized=False, fmpindex=fmp_index)

        else:
            eli_prod = lambda delta_wk: eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk, linearized=False, fmpindex=fmp_index)

    elif product == "SUM":
        if(hasDLRMesh): 
            raise NotImplementedError(
                "There is no implementation of the eliashberg product "
                "called %s when using DLR. Please use the FFT product instead."%product
            )

        eli_prod = lambda delta_wk: eliashberg_product(Gamma_pp_wk, g_wk, delta_wk, linearized=False, fmpindex=fmp_index)

    else:
        raise NotImplementedError(
            "There is no implementation of the eliashberg product"
            " called %s." % product
        )

    if not initial_delta:
        initial_delta = semi_random_initial_delta(g_wk)

    delta_prev = initial_delta

    it = 1 
    while True:

        if steffensens_method:
            fDelta = delta_prev - eli_prod(delta_prev)
            
            nonzeroinds = ~np.isclose(fDelta.data[:],0.0)
            gDelta = np.ones_like(fDelta.data[:])
            gDelta[nonzeroinds] = (delta_prev.data[nonzeroinds] + fDelta.data[nonzeroinds] - eli_prod(delta_prev + fDelta).data[nonzeroinds]) / fDelta.data[nonzeroinds] - 1.0

            delta_wk = delta_prev - fDelta / gDelta
        else:
            delta_wk = eli_prod(delta_prev)
        
        # Mixing of previous solution to avoid trapping
        delta_wk = mixing_frac * delta_prev + (1.0 - mixing_frac) * delta_wk

        # Symmetrizing
        delta_wk = symmetrize_fct(delta_wk)

        diff = np.max(np.abs(delta_prev.data[:] - delta_wk.data[:]))
        if verbose and (it%10==0):
            print("(%i) Error: %.14f"%(it, diff))
        if np.allclose(diff, 0, atol=tol):
            break

        delta_prev = delta_wk
        it += 1
        if it > max_it:
            raise AssertionError("Did not converge.")

    return delta_wk
