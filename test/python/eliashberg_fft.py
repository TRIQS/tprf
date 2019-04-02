# ----------------------------------------------------------------------

""" Compare final result and between steps of the Eliashberg equation to
previous erstablished results """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from pytriqs.gf import Gf, MeshImFreq, Idx

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lattice_dyson_g0_wk, solve_rpa_PH
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.lattice import gamma_PP_singlet
from triqs_tprf.lattice import eliashberg_product
from triqs_tprf.lattice import eliashberg_product_fft
from triqs_tprf.lattice import eliashberg_g_delta_g_product
from triqs_tprf.eliashberg import solve_eliashberg

# ----------------------------------------------------------------------

from utilities import read_TarGZ_HDFArchive

# ----------------------------------------------------------------------

# -- Read previous established results

import os
os.system('pwd')

filename = './eliashberg_benchmark.tar.gz'
p = read_TarGZ_HDFArchive(filename)['p']

print('\nThe benchmark data was obtained with the version of hash %s.'%p.tprf_hash)


t = - 2.0 * np.eye(1)
t_r = TBLattice(
    units = [(1, 0, 0)],
    hopping = {
        (+1,): t,
        (-1,): t,
        },
    orbital_positions = [(0,0,0)],
    orbital_names = ['up', 'do'],
    )

e_k = t_r.on_mesh_brillouin_zone((4, 1, 1))
print e_k.data
print p.e_k.data

#np.testing.assert_allclose(e_k.data, p.e_k.data)

# -- Test the construction of Gamma

gamma = gamma_PP_singlet(p.chi_c, p.chi_s, p.U_c, p.U_s)
np.testing.assert_allclose(gamma.data, p.gamma.data)

# -- Construct a Gamma with a possibly different wmesh than the benchmark data

print p.nw_gf
print p.nw_chi0

factor = 1.

wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(factor*p.nw_gf))
#g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=p.e_k, mesh=wmesh)
g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

wmesh_big = MeshImFreq(beta=p.beta, S='Fermion', n_max=int(2*factor*p.nw_gf))
g0_wk_big = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh_big)

#np.testing.assert_allclose(g0_wk.data, p.g0_wk.data)

#chi00_wk = imtime_bubble_chi0_wk(g0_wk_big, nw=int(factor*p.nw_chi0))
chi00_wk = imtime_bubble_chi0_wk(g0_wk_big, nw=int(2*factor*p.nw_gf + 1))
#chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=int(2*factor*p.nw_gf + 1))

print g0_wk.data.shape
print chi00_wk.data.shape
#exit()

chi_s = solve_rpa_PH(chi00_wk, +p.U_s)
chi_c = solve_rpa_PH(chi00_wk, -p.U_c) # Minus for correct charge rpa equation

gamma = gamma_PP_singlet(chi_c, chi_s, p.U_c, p.U_s)

# -- Test the Eliashberg equation with the new gamma

F_wk = eliashberg_g_delta_g_product(g0_wk, g0_wk)
next_delta = eliashberg_product(gamma, g0_wk, g0_wk)
next_delta_fft = eliashberg_product_fft(gamma, g0_wk, g0_wk)

from pytriqs.plot.mpl_interface import oplot, plt
subp = [3, 3, 1]
plt.figure(figsize=(9, 6))

plt.subplot(*subp); subp[-1] += 1
oplot(g0_wk[:, Idx(0, 0, 0)], label='d')
plt.title('g0_wk')

plt.subplot(*subp); subp[-1] += 1
oplot(chi00_wk[:, Idx(0, 0, 0)], label='d')
plt.title('chi00_wk')

plt.subplot(*subp); subp[-1] += 1
oplot(gamma[:, Idx(0, 0, 0)], label='d')
plt.title('gamma')

plt.subplot(*subp); subp[-1] += 1
oplot(F_wk[:, Idx(0, 0, 0)])
plt.title('F_wk')

plt.subplot(*subp); subp[-1] += 1
oplot(next_delta[:, Idx(0, 0, 0)], label='d')
plt.title('next_delta')

plt.subplot(*subp); subp[-1] += 1
oplot(next_delta_fft[:, Idx(0, 0, 0)], label='dfft')
plt.title('next_delta_fft')

plt.tight_layout()
plt.show(); exit()


print next_delta.data[:4, 0, 0, 0]
print next_delta_fft.data[:4, 0, 0, 0]

exit()

Es, eigen_modes = solve_eliashberg(gamma, g0_wk)

np.testing.assert_allclose(next_delta.data, p.next_delta.data, atol=1e-7)

np.testing.assert_allclose(Es[0], p.E)
try:
    np.testing.assert_allclose(eigen_modes[0].data, p.eigen_mode.data, atol=1e-6)
except AssertionError:
    np.testing.assert_allclose(-eigen_modes[0].data, p.eigen_mode.data, atol=1e-6)

print('\nSame results for both versions.')
