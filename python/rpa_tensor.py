
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.OperatorUtils import quartic_tensor_from_operator
from pyed.OperatorUtils import operator_from_quartic_tensor

from pyed.OperatorUtils import symmetrize_quartic_tensor

from pyed.OperatorUtils import quartic_permutation_symmetrize
from pyed.OperatorUtils import quartic_conjugation_symmetrize
from pyed.OperatorUtils import quartic_pauli_symmetrize

# ----------------------------------------------------------------------
def fundamental_operators_from_gf_struct(gf_struct):

    fundamental_operators = []
    for block, idxs in gf_struct:
        for idx in idxs:
            fundamental_operators.append( c(block, idx) )

    return fundamental_operators

# ----------------------------------------------------------------------
def get_rpa_tensor(H_int, fundamental_operators):

    U_int = quartic_tensor_from_operator(H_int, fundamental_operators)
    U_int = 4 * quartic_permutation_symmetrize(U_int)

    # -- Group in Gamma order cc^+cc^+ ( from c^+c^+cc )
    
    U_abcd = np.zeros_like(U_int)
    for a, b, c, d in itertools.product(range(U_int.shape[0]), repeat=4):
        U_abcd[a, b, c, d] = U_int[b, d, a, c]

    return U_abcd

# ----------------------------------------------------------------------
def get_gamma_rpa(chi0_wnn, U_abcd):

    # -- Build constant gamma

    gamma_wnn = chi0_wnn.copy()
    
    # Nb! In the three frequency form $\Gamma \propto U/\beta^2$

    beta = chi0_wnn.mesh.components[0].beta
    gamma_wnn.data[:] = U_abcd[None, None, None, ...] / beta**2
    
    return gamma_wnn

# ----------------------------------------------------------------------
def get_rpa_us_tensor(norbs, U, Up, J, Jp):
    """Returns interaction tensor for spin susceptibility 
    c^+cc^+c order

    Parameters:

    norbs : int, number of orbitals of the model
    U : float, intra orbital interaction
    Up : float, inter orbital interaction
    J : float, Hund's coupling
    Jp : float, pair hopping and spin flip
    ...
    """
    U_abcd = np.zeros(4*[norbs], dtype=complex)
    for a,b,c,d in itertools.product(range(norbs), repeat=4):
        if a == b == c == d:
            U_abcd[a,b,c,d] = U
        elif a == c != b == d:
            U_abcd[a,b,c,d] = Up
        elif a == b != c == d:
            U_abcd[a,b,c,d] = J
        elif a == d != b == c:
            U_abcd[a,b,c,d] = Jp

    return U_abcd
    
# ----------------------------------------------------------------------
def get_rpa_uc_tensor(norbs, U, Up, J, Jp):
    """Returns interaction tensor for charge susceptibility 
    c^+cc^+c order

    Parameters:

    norbs : int, number of orbitals of the model
    U : float, intra orbital interaction
    Up : float, inter orbital interaction
    J : float, Hund's coupling
    Jp : float, pair hopping and spin flip
    ...
    """
    U_abcd = np.zeros(4*[norbs], dtype=complex)
    for a,b,c,d in itertools.product(range(norbs), repeat=4):
        if a == b == c == d:
            U_abcd[a,b,c,d] = U
        elif a == c != b == d:
            U_abcd[a,b,c,d] = -Up + 2*J
        elif a == b != c == d:
            U_abcd[a,b,c,d] = 2*Up - J
        elif a == d != b == c:
            U_abcd[a,b,c,d] = Jp

    return U_abcd

# ----------------------------------------------------------------------    
def kanamori_charge_and_spin_quartic_interaction_tensors(norb, U, Up, J, Jp):

    """ Following Eliashberg notes, but with c+c+ cc order. """

    shape = [norb]*4
    U_c, U_s = np.zeros(shape), np.zeros(shape)
    
    for a, abar, b, bbar in itertools.product(range(norb), repeat=4):

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
            
        U_c[a, b, bbar, abar] = scalar_c
        U_s[a, b, bbar, abar] = scalar_s

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
    shape_8 = np.vstack(([2]*4, shape_4 / 2)).T.flatten()

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

    np.testing.assert_array_almost_equal(U_8[0, 0, 1, 1], zeros)
    np.testing.assert_array_almost_equal(U_8[1, 1, 0, 0], zeros)

    # -- split in charge and spin
    
    # c+ c c+ c form of the charge, spin diagonalization
    #U_c = U_8[0, 0, 0, 0] + U_8[0, 0, 1, 1]
    #U_s = U_8[0, 0, 0, 0] - U_8[0, 0, 1, 1]

    # c+ c+ c c  form of the charge, spin diagonalization
    U_c = U_8[0, 0, 0, 0] + U_8[0, 1, 1, 0]
    U_s = -U_8[0, 0, 0, 0] + U_8[0, 1, 1, 0]

    U_c *= 4
    U_s *= 4
    
    return U_c, U_s

# ----------------------------------------------------------------------    
def quartic_tensor_from_charge_and_spin(U_c, U_s):

    shape_4 = U_c.shape
    shape_8 = [2]*4 + list(shape_4)

    U_8 = np.zeros(shape_8, dtype=U_c.dtype)

    U_uu = 0.5 * (U_c - U_s)
    U_ud = 0.5 * (U_c + U_s)
    
    for s1, s2 in itertools.product(range(2), repeat=2):
        if s1 == s2:
            U_8[s1, s1, s1, s1] = U_uu
        else:
            U_8[s1, s2, s2, s1] = U_ud
            U_8[s1, s2, s1, s2] = -U_s

    U_8 = np.transpose(U_8, (0, 4, 1, 5, 2, 6, 3, 7))
    U_4 = U_8.reshape(2*np.array(shape_4))
    
    return 0.25*U_4
