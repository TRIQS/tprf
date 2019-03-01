# ----------------------------------------------------------------------

""" Comparison of the general RPA formalism and the matrix RPA formalism
for the spin- and charge-susceptibility. """

# ----------------------------------------------------------------------

import numpy as np

# ======================================================================
# General RPA formalism

# ----------------------------------------------------------------------
# -- One dimensional lattice with nearest neighbour hopping with
# -- explicit spin indice

t = 1.0
    
norb = 4
h_loc = np.zeros((norb, norb))
T = - t * np.eye(norb)

from triqs_tprf.tight_binding import TBLattice

t_r = TBLattice(
    units = [(1, 0, 0)],
    hopping = {
        # nearest neighbour hopping -t
        (0,): h_loc,
        (+1,): T,
        (-1,): T,
        },
    orbital_positions = [(0,0,0)] * norb,
    orbital_names = ['up_0', 'do_0', 'up_1', 'do_1'],
    )

e_k = t_r.on_mesh_brillouin_zone(n_k=(256, 1, 1))

# ----------------------------------------------------------------------
# -- Bare susceptibility from Green's function bubble

nw = 20
beta = 0.544

from pytriqs.gf import MeshImFreq, Idx
wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)

from triqs_tprf.lattice import lattice_dyson_g0_wk
g0_wk = lattice_dyson_g0_wk(mu=0., e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)
print
print 'chi0_q0 =\n', chi00_wk[Idx(0), Idx(0, 0, 0)].real.reshape((16,16))

# ----------------------------------------------------------------------
# -- Kanamori interaction

U = 2.4
J = 0.6

from pytriqs.operators.util.U_matrix import U_matrix_kanamori, U_matrix
from pytriqs.operators.util.hamiltonians import h_int_kanamori

orb_names = [0, 1]
spin_names = ['up', 'do']

U_ab, UPrime_ab = U_matrix_kanamori(
    n_orb=len(orb_names), U_int=U, J_hund=J)
    
H_int = h_int_kanamori(
    spin_names, orb_names, U_ab, UPrime_ab, J_hund=J,
    off_diag=False, map_operator_structure=None, H_dump=None) # orbital diag

print 'H_int =\n', H_int

# ----------------------------------------------------------------------
# -- RPA rank 4 (antisymmetrized) interaction tensor

from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.rpa_tensor import fundamental_operators_from_gf_struct

gf_struct = [['up_0', [0]], ['do_0',[0]], ['up_1', [0]], ['do_1', [0]]]

fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)
U_abcd = get_rpa_tensor(H_int, fundamental_operators)

print 'U_abdc =\n', U_abcd.reshape((16, 16)).real

# ----------------------------------------------------------------------
# -- Lattice susceptbility in the RPA approximation

from triqs_tprf.lattice import solve_rpa_PH
chi_wk = solve_rpa_PH(chi00_wk, U_abcd)

print 'chi_q0 =\n', chi_wk[Idx(0), Idx(0, 0, 0)].real.reshape((16,16))

Sz = 0.5 * np.diag([+1., -1., +1., -1.])
chi_SzSz_wk = chi_wk[0,0,0,0].copy()
chi_SzSz_wk.data[:] = np.einsum('wkabcd,ab,cd->wk', chi_wk.data, Sz, Sz)

n = np.diag([1., 1., 1., 1.])
chi_nn_wk = chi_wk[0,0,0,0].copy()
chi_nn_wk.data[:] = np.einsum('wkabcd,ab,cd->wk', chi_wk.data, n, n)

# ======================================================================
# Matrix RPA formalism

# ----------------------------------------------------------------------
# -- One dimensional lattice with nearest neighbour hopping without
# -- spin indice
    
norb = 2
h_loc = np.zeros((norb, norb))
T = - t * np.eye(norb)

t_r = TBLattice(
    units = [(1, 0, 0)],
    hopping = {
        # nearest neighbour hopping -t
        (0,): h_loc,
        (+1,): T,
        (-1,): T,
        },
    orbital_positions = [(0,0,0)] * norb,
    orbital_names = ['0', '1'],
    )

e_k = t_r.on_mesh_brillouin_zone(n_k=(256, 1, 1))

# ----------------------------------------------------------------------
# -- Bare susceptibility from Green's function bubble

from pytriqs.gf import MeshImFreq, Idx
wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)

from triqs_tprf.lattice import lattice_dyson_g0_wk
g0_wk = lattice_dyson_g0_wk(mu=0., e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)

print
print 'chi0_q0 =\n', chi00_wk[Idx(0), Idx(0, 0, 0)].real.reshape((4,4))

# ----------------------------------------------------------------------
# -- Spin- and charge-susceptibility in matrix RPA formalism

from triqs_tprf.lattice import solve_rpa_spin, solve_rpa_charge
from triqs_tprf.rpa_tensor import get_rpa_us_tensor, get_rpa_uc_tensor

Up = U - 2*J
Jp = J

U_abcd_spin = get_rpa_us_tensor(norb, U, Up, J, Jp)
U_abcd_charge = get_rpa_uc_tensor(norb, U, Up, J, Jp)

chi_spin_wk = solve_rpa_spin(chi00_wk, U_abcd_spin)
chi_charge_wk = solve_rpa_charge(chi00_wk, U_abcd_charge)

chi_SzSz_wk_ref = chi_spin_wk[0,0,0,0].copy()
chi_SzSz_wk_ref.data[:] = np.einsum('wkaabb->wk', chi_spin_wk.data)

chi_nn_wk_ref = chi_charge_wk[0,0,0,0].copy()
chi_nn_wk_ref.data[:] = np.einsum('wkaabb->wk', chi_charge_wk.data)

# ----------------------------------------------------------------------
# -- Compare the results of both formalisms. Note that in the matrix RPA
# -- the susceptibilties carry a factor of 0.5, 
# -- i.e. \chi_\text{general} = 2 * \chi_\text{matrix}.
# -- But because we build \chi^{S_zS_z} with the proper spin operator,
# -- which has \pm 0.5 elements we get an addtionl 0.25 factor.
# -- Therefore:
# -- \chi_\text{general}^{S_zS_z} = 0.5 * \chi_\text{matrix}^{S_zS_z},
# -- \chi_\text{general}^{nn} = 2 * \chi_\text{matrix}^{nn},

np.testing.assert_almost_equal(chi_SzSz_wk.data, 0.5*chi_SzSz_wk_ref.data)
np.testing.assert_almost_equal(chi_nn_wk.data, 2*chi_nn_wk_ref.data)
