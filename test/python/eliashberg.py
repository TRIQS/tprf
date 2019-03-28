# ----------------------------------------------------------------------

""" Compare final result and between steps of the Eliashberg equation to
previous erstablished results """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, MeshImFreq

from triqs_tprf.lattice import lattice_dyson_g0_wk, solve_rpa_PH
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.lattice import gamma_PP_singlet, eliashberg_product
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

from utilities import read_TarGZ_HDFArchive

# ----------------------------------------------------------------------

# -- Read previous established results

filename = './eliashberg_benchmark.tar.gz'
p = read_TarGZ_HDFArchive(filename)['p']

print('\nThe benchmark data was obtained with the version of hash %s.'%p.tprf_hash)

# -- Test the construction of Gamma

gamma = gamma_PP_singlet(p.chi_c, p.chi_s, p.U_c, p.U_s)
np.testing.assert_allclose(gamma.data, p.gamma.data)

# -- Construct a Gamma with a possibly different wmesh than the benchmark data

factor = 2.5

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(factor*p.nw_gf))
g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=p.e_k, mesh=wmesh)


chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=int(factor*p.nw_chi0))

chi_s = solve_rpa_PH(chi00_wk, p.U_s)
chi_c = solve_rpa_PH(chi00_wk, -p.U_c) # Minus for correct charge rpa equation

gamma = gamma_PP_singlet(chi_c, chi_s, p.U_c, p.U_s)

# -- Test the Eliashberg equation with the new gamma

next_delta = eliashberg_product(gamma, p.g0_wk, p.g0_wk) 
Es, eigen_modes = solve_eliashberg(gamma, p.g0_wk)

np.testing.assert_allclose(next_delta.data, p.next_delta.data, atol=10**(-7))

np.testing.assert_allclose(Es[0], p.E)
try:
    np.testing.assert_allclose(eigen_modes[0].data, p.eigen_mode.data, atol=10**(-6))
except AssertionError:
    np.testing.assert_allclose(-eigen_modes[0].data, p.eigen_mode.data, atol=10**(-6))

print('\nSame results for both versions.')
