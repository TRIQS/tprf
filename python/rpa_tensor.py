
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

