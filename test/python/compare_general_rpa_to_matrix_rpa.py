# ----------------------------------------------------------------------

""" Comparison of the general RPA formalism and the matrix RPA formalism
for the spin- and charge-susceptibility. """

# ----------------------------------------------------------------------

import numpy as np

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

from pytriqs.gf import MeshImFreq
from triqs_tprf.lattice import lattice_dyson_g0_wk

wmesh = MeshImFreq(beta=5.0, S='Fermion', n_max=30)
g0_wk = lattice_dyson_g0_wk(mu=0., e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)

# ----------------------------------------------------------------------
# -- Kanamori interaction

from pytriqs.operators.util import U_matrix_kanamori, h_int_kanamori
from pyed.OperatorUtils import fundamental_operators_from_gf_struct
from triqs_tprf.rpa_tensor import get_rpa_tensor
# -- Kanamori interaction

U = 1.0
J = 0.1

spin_names = ['up', 'do']
orb_names = range(norb)

# TRIQS uses spin as slow index
gf_struct = [ [spin_name, orb_names] for spin_name in spin_names ]
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

from triqs_tprf.matrix_rpa import lose_spin_degree_of_freedom, tprf_order_to_matrix_rpa_order

chi00_wk_wo_spin = lose_spin_degree_of_freedom(chi00_wk.data, rank=4, spin_fast=False)
chi00_wk_matrix_rpa = tprf_order_to_matrix_rpa_order(chi00_wk_wo_spin)

from triqs_tprf.matrix_rpa import get_rpa_us_tensor, get_rpa_uc_tensor

U = 1.0
Up = 0.8
J = 0.1
Jp = 0.1

us = get_rpa_us_tensor(norb, U, Up, J ,Jp)
uc = get_rpa_uc_tensor(norb, U, Up, J ,Jp)

from triqs_tprf.matrix_rpa import chi_rpa_spin, chi_rpa_charge

chi_spin = chi_rpa_spin(chi00_wk_matrix_rpa, us)
chi_charge = chi_rpa_charge(chi00_wk_matrix_rpa, uc)

uu = tprf_order_to_matrix_rpa_order(chi_wk.data[:,:,:norb,:norb,:norb,:norb])
ud = tprf_order_to_matrix_rpa_order(chi_wk.data[:,:,:norb,:norb,norb:,norb:])

np.testing.assert_allclose(uu, 0.5*(chi_charge + chi_spin))
np.testing.assert_allclose(ud, 0.5*(chi_charge - chi_spin))

np.testing.assert_allclose(chi_spin, uu-ud)
np.testing.assert_allclose(chi_charge, uu+ud)

from triqs_tprf.rpa_tensor import split_quartic_tensor_in_charge_and_spin
uc_ref, us_ref = split_quartic_tensor_in_charge_and_spin(U_abcd)

np.testing.assert_allclose(us, tprf_order_to_matrix_rpa_order(us_ref))
np.testing.assert_allclose(uc, tprf_order_to_matrix_rpa_order(uc_ref))
