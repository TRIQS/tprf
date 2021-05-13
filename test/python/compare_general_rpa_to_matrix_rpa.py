# ----------------------------------------------------------------------

""" Comparison of the general RPA formalism and the matrix RPA formalism
for the spin- and charge-susceptibility. 

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import numpy as np
import itertools

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

# ======================================================================
# General RPA formalism

# ----------------------------------------------------------------------
# -- Two dimensional lattice with nearest neighbour hopping intra and
# -- inter orbital and explicit spin indice

norb = 2
spin_names = ['up', 'do']

n_orb_spin = norb * len(spin_names)

t_intra = 1.0
t_inter = 0.1

inter_orbital_hopping = np.zeros((n_orb_spin, n_orb_spin))
inter_orbital_hopping[0,1] = inter_orbital_hopping[1,0] = 1
inter_orbital_hopping[2,3] = inter_orbital_hopping[3,2] = 1

H = TBLattice(
    units = [(1, 0, 0), (0, 1, 0)],
    hopping = {
        # nearest neighbour hopping -t
        ( 0,+1): -t_intra * np.eye(n_orb_spin) - t_inter * inter_orbital_hopping,
        ( 0,-1): -t_intra * np.eye(n_orb_spin) - t_inter * inter_orbital_hopping,
        (+1, 0): -t_intra * np.eye(n_orb_spin) - t_inter * inter_orbital_hopping,
        (-1, 0): -t_intra * np.eye(n_orb_spin) - t_inter * inter_orbital_hopping,
        },
    orbital_positions = [(0,0,0)]*n_orb_spin,
    )

e_k = H.on_mesh_brillouin_zone(n_k = (16, 16, 1))
print(e_k)

# ----------------------------------------------------------------------
# -- Bare susceptibility from Green's function bubble

from triqs.gf import MeshImFreq
from triqs_tprf.lattice import lattice_dyson_g0_wk

wmesh = MeshImFreq(beta=5.0, S='Fermion', n_max=30)
g0_wk = lattice_dyson_g0_wk(mu=0., e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)

# ----------------------------------------------------------------------
# -- Kanamori interaction

from triqs.operators.util import U_matrix_kanamori, h_int_kanamori
from triqs_tprf.OperatorUtils import fundamental_operators_from_gf_struct
from triqs_tprf.rpa_tensor import get_rpa_tensor

U = 1.0
J = 0.1

spin_names = ['up', 'do']
orb_names = list(range(norb))

# TRIQS uses spin as slow index
gf_struct = [ [spin_name, norb] for spin_name in spin_names ]
Umat, Upmat = U_matrix_kanamori(n_orb=norb, U_int=U, J_hund=J)
H_int = h_int_kanamori(spin_names, orb_names, U=Umat, Uprime=Upmat, J_hund=J, off_diag=True)

fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)

U_abcd = get_rpa_tensor(H_int, fundamental_operators) # given in cc^+cc^+

# ----------------------------------------------------------------------
# -- Lattice susceptbility in the RPA approximation

from triqs_tprf.lattice import solve_rpa_PH

chi_wk = solve_rpa_PH(chi00_wk, U_abcd)

# ======================================================================
# Matrix RPA formalism

# ----------------------------------------------------------------------
# -- Showcase reshaping from 4-rank tensors to matrix as done in papers
# -- and test if the process is reversable

from triqs_tprf.matrix_rpa import tensor_to_matrix, matrix_to_tensor

test_chi = np.chararray([norb]*4, itemsize=4)

for i,j,k,l in itertools.product(list(range(norb)), repeat=4):
    test_chi[i,j,k,l] = "%s%s%s%s"%(i,j,k,l)

print()
print(tensor_to_matrix(test_chi))
print()
np.testing.assert_equal(test_chi, matrix_to_tensor(tensor_to_matrix(test_chi)))

test_us = np.chararray([norb]*4, itemsize=2)

for a,b,c,d in itertools.product(list(range(norb)), repeat=4):

    if a == b == c == d:
        test_us[a,b,c,d] = 'U'

    elif a == c != b == d:
        test_us[a,b,c,d] = 'Up'

    elif a == b != c == d:
        test_us[a,b,c,d] = 'J'

    elif a == d != b == c:
        test_us[a,b,c,d] = 'Jp'

    else:
        test_us[a,b,c,d] = '0'
print()
print(tensor_to_matrix(test_us))
print()
np.testing.assert_equal(test_us, matrix_to_tensor(tensor_to_matrix(test_us)))

test_uc = np.chararray([norb]*4, itemsize=6)

for a,b,c,d in itertools.product(list(range(norb)), repeat=4):

    if a == b == c == d:
        test_uc[a,b,c,d] = 'U'

    elif a == c != b == d:
        test_uc[a,b,c,d] = '-Up+2J'

    elif a == b != c == d:
        test_uc[a,b,c,d] = '2Up-J'

    elif a == d != b == c:
        test_uc[a,b,c,d] = 'Jp'

    else:
        test_uc[a,b,c,d] = '0'

print()
print(tensor_to_matrix(test_uc))
print()
np.testing.assert_equal(test_uc, matrix_to_tensor(tensor_to_matrix(test_uc)))

# ----------------------------------------------------------------------
# -- Calculate chi spin/charge with spin independent chi0

from triqs_tprf.matrix_rpa import lose_spin_degree_of_freedom, tprf_order_to_matrix_rpa_order

chi00_wk_wo_spin_array = lose_spin_degree_of_freedom(chi00_wk.data, 
                                                        rank=4, spin_fast=False) # c^+cc^+c
chi00_wk_matrix_rpa = tprf_order_to_matrix_rpa_order(chi00_wk_wo_spin_array) # now in c^+ccc^+ order

from triqs_tprf.matrix_rpa import get_rpa_us_tensor, get_rpa_uc_tensor

U = 1.0
Up = 0.8
J = 0.1
Jp = 0.1

us_matrix_rpa = get_rpa_us_tensor(norb, U, Up, J ,Jp) # given in cc^+c^+c
uc_matrix_rpa = get_rpa_uc_tensor(norb, U, Up, J ,Jp) # given in cc^+c^+c

from triqs_tprf.matrix_rpa import chi_rpa_spin, chi_rpa_charge

chi_spin_matrix_rpa = chi_rpa_spin(chi00_wk_matrix_rpa, us_matrix_rpa) # given in c^+ccc^+
chi_charge_matrix_rpa = chi_rpa_charge(chi00_wk_matrix_rpa, uc_matrix_rpa) # given in c^+ccc^+

# ======================================================================
# Test if the TPRF and matrix RPA calculations are compatible 

# -- Construct the susceptibiilites with spin up up and up down
chi_uu = tprf_order_to_matrix_rpa_order(chi_wk.data[:,:,:norb,:norb,:norb,:norb])
chi_ud = tprf_order_to_matrix_rpa_order(chi_wk.data[:,:,:norb,:norb,norb:,norb:])

# -- Test Equation \chi_{\uparrow \uparrow} = 0.5 * (\chi^{c} + \chi^{s})
np.testing.assert_allclose(chi_uu, 0.5*(chi_charge_matrix_rpa + chi_spin_matrix_rpa))
# -- Test Equation \chi_{\uparrow \downarrow} = 0.5 * (\chi^{c} - \chi^{s})
np.testing.assert_allclose(chi_ud, 0.5*(chi_charge_matrix_rpa - chi_spin_matrix_rpa))

# -- Test Equation \chi^{s} = \chi_{\uparrow \uparrow} - \chi_{\uparrow \downarrow}
np.testing.assert_allclose(chi_spin_matrix_rpa, chi_uu-chi_ud)
# -- Test Equation \chi^{c} = \chi_{\uparrow \uparrow} + \chi_{\uparrow \downarrow}
np.testing.assert_allclose(chi_charge_matrix_rpa, chi_uu+chi_ud)

# -- Construct the susceptibiilites with spin up down but crossed
chi_xud = tprf_order_to_matrix_rpa_order(chi_wk.data[:,:,:norb,norb:,norb:,:norb])
# -- Test Equation \chi_{\times \uparrow \downarrow} = \chi^{s}
np.testing.assert_allclose(chi_spin_matrix_rpa, chi_xud)

# -- Obtain the charge and spin interaction tensor via the splitting of the quartic tensor
from triqs_tprf.rpa_tensor import split_quartic_tensor_in_charge_and_spin
uc, us = split_quartic_tensor_in_charge_and_spin(U_abcd) # Given in cc^+cc^+

# -- Order cc^+cc^+ -> cc^+c^+c and compare
np.testing.assert_allclose(us_matrix_rpa, tprf_order_to_matrix_rpa_order(us))
np.testing.assert_allclose(uc_matrix_rpa, tprf_order_to_matrix_rpa_order(uc))

# ======================================================================
# Use spin-independent TPRF formalism

from triqs_tprf.rpa_tensor import lose_spin_degree_of_freedom

chi00_wk_wo_spin = lose_spin_degree_of_freedom(chi00_wk, spin_fast=False)

from triqs_tprf.rpa_tensor import general_susceptibility_from_charge_and_spin

chi_spin = solve_rpa_PH(chi00_wk_wo_spin, us)
chi_charge = solve_rpa_PH(chi00_wk_wo_spin, -uc)

chi_wk_ref = general_susceptibility_from_charge_and_spin(chi_charge, chi_spin, spin_fast=False)

# -- Test if spin/charge susceptibilites are the same with matrix RPA
np.testing.assert_allclose(chi_spin_matrix_rpa, tprf_order_to_matrix_rpa_order(chi_spin.data))
np.testing.assert_allclose(chi_charge_matrix_rpa, tprf_order_to_matrix_rpa_order(chi_charge.data))

# -- Test if general susceptibility is unchanged for a spin-independent calculation
np.testing.assert_allclose(chi_wk.data, chi_wk_ref.data)

# -- Test if general -> spin/charge -> general is the same

from triqs_tprf.rpa_tensor import charge_and_spin_susceptibility_from_general

chi_charge_ref, chi_spin_ref = charge_and_spin_susceptibility_from_general(chi_wk, spin_fast=False)
chi_wk_ref = general_susceptibility_from_charge_and_spin(chi_charge_ref, chi_spin_ref,
                                                                            spin_fast=False)

np.testing.assert_allclose(chi_wk.data, chi_wk_ref.data)
