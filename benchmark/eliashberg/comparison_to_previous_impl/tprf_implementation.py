# ----------------------------------------------------------------------

""" Compare the TPRF implementation of the Eliashberg equation to
a previous implementation, which results are stored in a file.
"""

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from h5 import HDFArchive
import triqs.utility.mpi as mpi

from triqs.lattice.lattice_tools import BravaisLattice as BravaisLattice
from triqs.lattice.lattice_tools import BrillouinZone as BrillouinZone
from triqs.gf import Gf, MeshBrillouinZone, MeshImFreq

# ----------------------------------------------------------------------

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

# ================== Benchmark data ===========================================

# -- Load data from previous implementation of the eliashberg equation

INPUT_FILENAME = './eliashberg_benchmark_k_mesh_2.h5'

load = HDFArchive(INPUT_FILENAME)

# -- Get parameters

parameters = load['parameters']

beta = parameters['beta']
n_max = parameters['iw_mesh']
mu = parameters['chemical_potential']
U = parameters['U']

nk = parameters['k_mesh'][0]
n_max = load['chi_c'].shape[0]/2 + 1

# -- Load previously calculated quantities and reshape to (wmesh, kmesh)

e_k_ref = load['e_k']
chi_s_ref = load['chi_s'].reshape(-1, nk**2)
chi_c_ref = load['chi_c'].reshape(-1, nk**2)
gamma_ref = load['gamma'].reshape(-1, nk**2)
init_delta_ref = load['init_delta'].reshape(-1, nk**2)
next_delta_ref = load['next_delta'].reshape(-1, nk**2)
final_delta_ref = load['final_delta'].reshape(-1, nk**2)
lamb_ref = load['lambda']
negative_final_delta_ref = load['negative_final_delta'].reshape(-1, nk**2)
negative_lamb_ref = load['negative_lambda']

# ================== TPRF calculations ========================================

# -- Create dispersion relation Green's function object

norbs = e_k_ref.shape[-1]

units = [(1,0,0), (0,1,0), (0,0,1)]
orbital_positions = [(0, 0, 0)]
bl = BravaisLattice(units, orbital_positions)
bz = BrillouinZone(bl)
periodization_matrix = nk * np.eye(3, dtype=np.int32)
periodization_matrix[2,2] = 1 
kmesh = MeshBrillouinZone(bz, periodization_matrix)

e_k = Gf(mesh=kmesh, target_shape=[norbs, norbs])

e_k.data[:] = e_k_ref.reshape(nk**2, norbs, norbs)

# --  Calculate bare Green's function

wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_max)
g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

# -- Calculate bare bubble

chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=n_max)

# -- Calculate chi spin and charge

U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(norbs, U, 0, 0, 0)

chi_s = solve_rpa_PH(chi00_wk, U_s)
chi_c = solve_rpa_PH(chi00_wk, -U_c) # Minus for correct charge rpa equation

# For the first few comparisons the benchmark data is very close to the TPRF results
atol = 10**(-8)

# Compare chi spin and charge to benchmark data
np.testing.assert_allclose(chi_s.data.reshape(-1, nk**2), chi_s_ref, atol=atol)
np.testing.assert_allclose(chi_c.data.reshape(-1, nk**2), chi_c_ref, atol=atol)

# -- Create Gamma (singlet particle-particle vertex)

gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)

# Compare Gamma to benchmark data
np.testing.assert_allclose(gamma.data.reshape(-1, nk**2), gamma_ref, atol=atol)

# -- Use the same inital delta as was used in the previous implementation

init_delta = g0_wk.copy()
init_delta.data[:] = init_delta_ref.reshape(-1, nk**2, 1, 1)

# Compare inital deltas
np.testing.assert_allclose(init_delta.data.reshape(-1, nk**2), init_delta_ref, atol=atol)

# -- Use the Eliashberg product once on the inital delta

next_delta = eliashberg_product(gamma, g0_wk, init_delta) 

# For the next steps the benchmark data differs more from the TPRF results. This is
# properly due to the Fourier transformation version of the eliashberg product
# in the benchmark data and the use of different eigenvalue search algorithms.
atol = 10**(-6)

# Compare the output to the benchmark data
np.testing.assert_allclose(next_delta.data.reshape(-1, nk**2), next_delta_ref, atol=atol)

# -- Solve the Eliashberg equation for the maximum eigenvalue and its eigenvector

E, eigen_modes = solve_eliashberg(gamma, g0_wk)

# Compare the eigenvalue to the benchmark data
np.testing.assert_allclose(E[0], lamb_ref)

# Helper function to compare eigenvectors with a phase
def rotate_to_real_plane(arr):

    if np.max(np.abs(arr.imag)) < 10**(-8):
        return arr

    angles = np.angle(arr[np.abs(arr) > 10**(-8)])
    
    # Check if eigenvector has a global phase
    np.testing.assert_allclose(np.max(angles) - np.min(angles), np.pi, atol=atol)

    real_arr = np.exp(-1j * np.min(angles)) * arr
    
    # Check if array is now only real
    np.testing.assert_allclose(np.max(np.abs(real_arr.imag)), 0, atol=atol)

    return real_arr

final_delta = rotate_to_real_plane(eigen_modes[0].data).reshape(-1, nk*2)
final_delta_ref = rotate_to_real_plane(final_delta_ref)

# Compare the eigenvector to the benchmark data. Test for a  multiplication with -1
try:
    np.testing.assert_allclose(final_delta, final_delta_ref, atol=atol)
except AssertionError:
    np.testing.assert_allclose(-final_delta, final_delta_ref, atol=atol)

