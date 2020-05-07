# ----------------------------------------------------------------------

""" Comparison of RPA magnetic response at q=0 in a two band 
Kanamori model from TPRF and from analytic expressions. """

# ----------------------------------------------------------------------

import numpy as np

from scipy.integrate import quad
from scipy.optimize import brentq

np.set_printoptions(precision=3, linewidth=1000)

# ----------------------------------------------------------------------
# -- Analytic expressions for q=0 response of one dim model w nn hopping.

def get_density_of_states(eps, t):
    return 1./np.sqrt(1. - (eps/(2*t))**2) / (2*t * np.pi)

def fermi_distribution(beta, eps):
    return 1./(np.exp(beta*eps) + 1.)

def fermi_distribution_derivative(beta, eps):
    return -beta/4. / np.cosh(0.5*beta*eps)**2

def chi0_q0_integral(t, beta):

    def integrand(eps):
        rho = get_density_of_states(eps, t)
        df = fermi_distribution_derivative(beta, eps)
        return -rho * df

    chi0_q0, err = quad(integrand, -2.*t, 2.*t)

    return chi0_q0        

def find_Uc(t, beta):

    def root_function(U):
        chi0_q0 = chi0_q0_integral(t, beta)
        return 1 - U * chi0_q0

    Uc = brentq(root_function, 0, 100.)
    return Uc

# ----------------------------------------------------------------------
# -- One dimensional lattice with nearest neighbour hopping

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

from triqs.gf import MeshImFreq, Idx
wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)

from triqs_tprf.lattice import lattice_dyson_g0_wk
g0_wk = lattice_dyson_g0_wk(mu=0., e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)
print()
print('chi0_q0 =\n', chi00_wk[Idx(0), Idx(0, 0, 0)].real.reshape((16,16)))

#print
#import itertools
#for idxs in itertools.product(range(2), repeat=4):
#    print idxs, chi00_wk[Idx(0), Idx(0, 0, 0)].real[idxs]
    
chi0_q0_ref = chi0_q0_integral(t, beta)

print()
print('chi0_q0     =', chi00_wk[Idx(0), Idx(0, 0, 0)][0,0,0,0].real)
print('chi0_q0_ref =', chi0_q0_ref)
print()

# ----------------------------------------------------------------------
# -- Kanamori interaction

U = 2.4
J = 0.6

from triqs.operators.util.U_matrix import U_matrix_kanamori, U_matrix
from triqs.operators.util.hamiltonians import h_int_kanamori

orb_names = [0, 1]
spin_names = ['up', 'do']

U_ab, UPrime_ab = U_matrix_kanamori(
    n_orb=len(orb_names), U_int=U, J_hund=J)
    
H_int = h_int_kanamori(
    spin_names, orb_names, U_ab, UPrime_ab, J_hund=J,
    off_diag=False, map_operator_structure=None, H_dump=None) # orbital diag

print('H_int =\n', H_int)

# ----------------------------------------------------------------------
# -- RPA rank 4 (antisymmetrized) interaction tensor

from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.rpa_tensor import fundamental_operators_from_gf_struct

gf_struct = [['up_0', [0]], ['do_0',[0]], ['up_1', [0]], ['do_1', [0]]]

fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)
U_abcd = get_rpa_tensor(H_int, fundamental_operators)

print('U_abdc =\n', U_abcd.reshape((16, 16)).real)

# ----------------------------------------------------------------------
# -- Lattice susceptbility in the RPA approximation

from triqs_tprf.lattice import solve_rpa_PH
chi_wk = solve_rpa_PH(chi00_wk, U_abcd)

print('chi_q0 =\n', chi_wk[Idx(0), Idx(0, 0, 0)].real.reshape((16,16)))

#print
#for idxs in itertools.product(range(4), repeat=4):
#    print idxs, chi_wk[Idx(0), Idx(0, 0, 0)][idxs].real

Sz = 0.5 * np.diag([+1., -1., +1., -1.])
chi_SzSz_wk = chi_wk[0,0,0,0].copy()
chi_SzSz_wk.data[:] = np.einsum('wkabcd,ab,cd->wk', chi_wk.data, Sz, Sz)

# Eq. 7.43 Fazekas (additional 0.5 factor)
Ueff = (U + (norb/2 - 1)*J)
chi_q0_ref = chi0_q0_ref / (1. -  Ueff * chi0_q0_ref)

print()
print('chi_SzSz_q0 =', chi_SzSz_wk[Idx(0), Idx(0, 0, 0)].real)
print('chi_q0_ref =', chi_q0_ref)
print()

np.testing.assert_array_almost_equal(
    chi_SzSz_wk[Idx(0), Idx(0, 0, 0)].real, chi_q0_ref)
