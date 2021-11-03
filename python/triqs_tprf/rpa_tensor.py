# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019, The Simons Foundation and S. Käser
# Author: H. U.R. Strand, S. Käser
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

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Gf
from triqs.gf.block_gf import fix_gf_struct_type
from triqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from triqs_tprf.OperatorUtils import quartic_tensor_from_operator
from triqs_tprf.OperatorUtils import operator_from_quartic_tensor

from triqs_tprf.OperatorUtils import symmetrize_quartic_tensor

from triqs_tprf.OperatorUtils import quartic_permutation_symmetrize
from triqs_tprf.OperatorUtils import quartic_conjugation_symmetrize
from triqs_tprf.OperatorUtils import quartic_pauli_symmetrize

# ----------------------------------------------------------------------
def fundamental_operators_from_gf_struct(gf_struct):

    gf_struct = fix_gf_struct_type(gf_struct)

    fundamental_operators = []
    for block, block_size in gf_struct:
        for idx in range(block_size):
            fundamental_operators.append( c(block, idx) )

    return fundamental_operators

# ----------------------------------------------------------------------
def get_rpa_tensor(H_int, fundamental_operators):

    """ Takes a TRIQS operator object and extracts the quartic terms
    and returns the corresponding antisymmetrized quartic tensor in
    vertex index order, i.e., cc+cc+. """
    
    U_abcd = quartic_tensor_from_operator(H_int, fundamental_operators)
    U_abcd = 4 * quartic_permutation_symmetrize(U_abcd)

    # -- Group in Gamma order cc^+cc^+ ( from c^+c^+cc )
    Gamma_RPA_abcd = np.ascontiguousarray(np.transpose(U_abcd, (2, 0, 3, 1)))
        
    return Gamma_RPA_abcd

# ----------------------------------------------------------------------
def get_gamma_rpa(chi0_wnn, U_abcd):

    # -- Build constant gamma

    gamma_wnn = chi0_wnn.copy()
    
    # Nb! In the three frequency form $\Gamma \propto U/\beta^2$

    beta = chi0_wnn.mesh.components[0].beta
    gamma_wnn.data[:] = U_abcd[None, None, None, ...] / beta**2
    
    return gamma_wnn

# ----------------------------------------------------------------------    
def kanamori_charge_and_spin_quartic_interaction_tensors(norb, U, Up, J, Jp):

    """ Following Eliashberg notes. """

    shape = [norb]*4
    U_c, U_s = np.zeros(shape, dtype=complex), np.zeros(shape, dtype=complex)
    
    for a, abar, b, bbar in itertools.product(list(range(norb)), repeat=4):

        scalar_c, scalar_s = 0, 0

        if (a == abar) and (a == b) and (b == bbar):
            scalar_c = U
            scalar_s = U
            
        if (a == bbar) and  (a != b) and (abar == b):
            scalar_c = -Up + 2*J
            scalar_s = Up

        if (a == abar) and  (a != b) and (b == bbar):
            scalar_c = 2*Up - J
            scalar_s = J

        if (a == b) and  (a != abar) and (abar == bbar):
            scalar_c = Jp
            scalar_s = Jp
            
        #U_c[a, b, bbar, abar] = scalar_c
        #U_s[a, b, bbar, abar] = scalar_s
        
        U_c[a, abar, b, bbar] = scalar_c
        U_s[a, abar, b, bbar] = scalar_s

    return U_c, U_s

# ----------------------------------------------------------------------    
def split_quartic_tensor_in_charge_and_spin(U_4):

    """Assuming spin is the slow index and orbital is the fast index 

    Using a rank 4 U_abcd tensor with composite (spin, orbital) indices
    as input, assuming c^+c^+ c c structure of the tensor

    Returns:

    U_c : Charge channel rank 4 interaction tensor
    U_s : Spin channel rank 4 interaction tensor """

    shape_4 = np.array(U_4.shape)
    shape_8 = np.vstack(([2]*4, shape_4 // 2)).T.flatten()

    U_8 = U_4.reshape(shape_8)

    U_8 = np.transpose(U_8, (0, 2, 4, 6, 1, 3, 5, 7)) # spin first

    # -- Check spin-conservation

    zeros = np.zeros_like(U_8[0, 0, 0, 0])
    
    np.testing.assert_array_almost_equal(U_8[0, 0, 0, 1], zeros)
    np.testing.assert_array_almost_equal(U_8[0, 0, 1, 0], zeros)
    np.testing.assert_array_almost_equal(U_8[0, 1, 0, 0], zeros)
    np.testing.assert_array_almost_equal(U_8[1, 0, 0, 0], zeros)

    np.testing.assert_array_almost_equal(U_8[1, 1, 1, 0], zeros)
    np.testing.assert_array_almost_equal(U_8[1, 1, 0, 1], zeros)
    np.testing.assert_array_almost_equal(U_8[1, 0, 1, 1], zeros)
    np.testing.assert_array_almost_equal(U_8[0, 1, 1, 1], zeros)

    np.testing.assert_array_almost_equal(U_8[1, 0, 1, 0], zeros)
    np.testing.assert_array_almost_equal(U_8[0, 1, 0, 1], zeros)

    # -- split in charge and spin
    
    # c+ c c+ c form of the charge, spin diagonalization
    U_c = - U_8[0, 0, 0, 0] - U_8[0, 0, 1, 1]
    U_s = U_8[0, 0, 0, 0] - U_8[0, 0, 1, 1]

    # c+ c+ c c  form of the charge, spin diagonalization
    #U_c = U_8[0, 0, 0, 0] + U_8[0, 1, 1, 0]
    #U_s = -U_8[0, 0, 0, 0] + U_8[0, 1, 1, 0]

    #U_c *= 4
    #U_s *= 4
    
    return U_c, U_s

# ----------------------------------------------------------------------    
def quartic_tensor_from_charge_and_spin(U_c, U_s):

    shape_4 = U_c.shape
    shape_8 = [2]*4 + list(shape_4)

    U_8 = np.zeros(shape_8, dtype=U_c.dtype)

    U_uu = -0.5 * (U_c - U_s)
    U_ud = +0.5 * (U_c + U_s)
    
    for s1, s2 in itertools.product(list(range(2)), repeat=2):
        if s1 == s2:
            U_8[s1, s1, s1, s1] = U_uu
        else:
            U_8[s1, s1, s2, s2] = -U_ud
            U_8[s1, s2, s2, s1] = U_s

    U_8 = np.transpose(U_8, (0, 4, 1, 5, 2, 6, 3, 7))
    U_4 = U_8.reshape(2*np.array(shape_4))
    
    return U_4

# ----------------------------------------------------------------------    
def kanamori_quartic_tensor(norb, U, Up, J, Jp):
    r"""Return Kanamori interaction as a quartic tensor

    .. math::

        \hat{U}_{\text { Kanamori }} = U \sum_{i} \hat{n}_{i, \uparrow} \hat{n}_{i, \downarrow}+
        \sum_{i>j, s, s^{\prime}}\left(U^{\prime}-J \delta_{\sigma, \sigma^{\prime}}\right)
        \hat{n}_{i, \sigma} \hat{n}_{j, \sigma^{\prime}} - 
        \\ J \sum_{i \neq j}\left(\hat{c}_{i, \downarrow}^{\dagger} \hat{c}_{j, \uparrow}^{\dagger}
        \hat{c}_{j, \downarrow} \hat{c}_{i, \uparrow}+\hat{c}_{j, \uparrow}^{\dagger}
        \hat{c}_{j, \downarrow}^{\dagger} \hat{c}_{i, \uparrow} \hat{c}_{i, \downarrow}+
        \mathrm{h.c.}\right)

    Parameters
    ----------
    norb : int,
           Number of orbitals excluding spin.
    U : complex,
        Strength of intra-orbital interaction.
    Up : complex,
         Strength of inter-orbital interaction.
    J : complex, 
        Strength of Hound's coupling.
    Jp : complex,
         Strength pair hopping and spin-flip.

    Returns
    -------
    U : np.ndarray,
        shape = (2*norb, 2*norb, 2*norb, 2*norb)
    """

    U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(norb, U, Up, J, Jp)
    U = quartic_tensor_from_charge_and_spin(U_c, U_s)

    return U

# ----------------------------------------------------------------------    
def lose_spin_degree_of_freedom(gf, spin_fast=True):
    """Only keep the up spin elements of a Greens function

    Parameters:

    gf: Greens function, the last rank dimenions must the orbitals
    spin_fast: bool, True if spin is the fast index, e.g.
                    xz up, xz down, xy up, xy down, yz up, yz down,
                or False if spin is the slow index, e.g.
                    xz up, xy up, yz up, xz down, xy down, yz down.
    """
    
    norb = gf.target_shape[-1] // 2

    if spin_fast:
        idx = gf.target_rank*(slice(None, None, 2),)

    else:
        idx = gf.target_rank*(slice(norb, None),)

    return gf[idx]

# ----------------------------------------------------------------------    
def general_susceptibility_from_charge_and_spin(chi_c, chi_s, spin_fast=True):
    """Construct a general susceptibility (spin dependent) from chi spin and charge

    Parameters:

    chi_c: Greens function, the charge susceptibility
    chi_s: Greens function, the spin susceptibility
    spin_fast: bool, True if spin is the fast index, e.g.
                        xz up, xz down, xy up, xy down, yz up, yz down,
                    or False if spin is the slow index, e.g.
                        xz up, xy up, yz up, xz down, xy down, yz down.
    """

    norb = chi_c.target_shape[-1]
    rank = chi_c.rank
    target_rank = chi_c.target_rank

    chi_general = Gf(mesh=chi_c.mesh, target_shape=target_rank*(2*norb,))
    
    chi_uu = 0.5 * (chi_c + chi_s)
    chi_ud = 0.5 * (chi_c - chi_s)
    chi_xud = chi_s

    idx_rank = rank * (slice(None),)

    if spin_fast:
        up = slice(None, None, 2)
        down = slice(1, None, 2)

    else:
        up = slice(norb)
        down = slice(norb, None)

    chi_general.data[idx_rank + (up, up, up, up)] = chi_uu.data
    chi_general.data[idx_rank + (down, down, down, down)] = chi_uu.data
    chi_general.data[idx_rank + (up, up, down, down)] = chi_ud.data
    chi_general.data[idx_rank + (down, down, up, up)] = chi_ud.data
    chi_general.data[idx_rank + (up, down, down, up)] = chi_xud.data
    chi_general.data[idx_rank + (down, up, up, down)] = chi_xud.data

    return chi_general

# ----------------------------------------------------------------------    
def charge_and_spin_susceptibility_from_general(chi, spin_fast=True, check_spin_conservation=True):
    r"""Construct a chi spin and charge from a generalized susceptibility

    Should only be used for a :math:`SU(2)` susceptibility.

    Parameters
    ----------
    chi : Gf,
          Generalized susceptibility :math:`\chi_{a,b,c,d}` where :math:`a,b,c,d` are
          combined indices of spin and orbital.
    spin_fast : bool, optional
                True if spin is the fast index, e.g.
                xz up, xz down, xy up, xy down, yz up, yz down.
                False if spin is the slow index, e.g.
                xz up, xy up, yz up, xz down, xy down, yz down.
    check_spin_conservation : bool, optional
                              True if the susceptibility should be checked for spin
                              conservation, False otherwise.           
    """

    norb = chi.target_shape[-1] // 2
    rank = chi.rank
    target_rank = chi.target_rank

    idx_rank = rank * (slice(None),)

    if spin_fast:
        up = slice(None, None, 2)
        down = slice(1, None, 2)

    else:
        up = slice(norb)
        down = slice(norb, None)

    if check_spin_conservation:
        np.testing.assert_allclose(chi[(up, up, up, down)].data, 0)
        np.testing.assert_allclose(chi[(up, up, down, up)].data, 0)
        np.testing.assert_allclose(chi[(up, down, up, up)].data, 0)
        np.testing.assert_allclose(chi[(down, up, up, up)].data, 0)

        np.testing.assert_allclose(chi[(down, down, down, up)].data, 0)
        np.testing.assert_allclose(chi[(down, down, up, down)].data, 0)
        np.testing.assert_allclose(chi[(down, up, down, down)].data, 0)
        np.testing.assert_allclose(chi[(up, down, down, down)].data, 0)

        np.testing.assert_allclose(chi[(up, down, up, down)].data, 0)
        np.testing.assert_allclose(chi[(down, up, down, up)].data, 0)

    chi_uu = chi[(up, up, up, up)]
    chi_ud = chi[(up, up, down, down)]

    chi_s = chi_uu - chi_ud
    chi_c = chi_uu + chi_ud

    return chi_c, chi_s



