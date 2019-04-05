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
from triqs_tprf.lattice import split_into_dynamic_and_constant
from triqs_tprf.eliashberg import solve_eliashberg
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors

# ----------------------------------------------------------------------

from utilities import read_TarGZ_HDFArchive

# ----------------------------------------------------------------------

# -- Read previous established results

import os
os.system('pwd')

nk = 16
nw = 200
beta = 5
mu = 0.0
U = 1.0


t = - 2.0 * np.eye(1)
t_r = TBLattice(
    units = [(1, 0, 0)],
    hopping = {
        (+1,): t,
        (-1,): t,
        },
    orbital_positions = [(0,0,0)],
    )

e_k = t_r.on_mesh_brillouin_zone((nk, 1, 1))

factor = 2.

wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=int(nw))

g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)

wmesh_big = MeshImFreq(beta=beta, S='Fermion', n_max=int(2*factor*nw))
g0_wk_big = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh_big)

chi00_wk = imtime_bubble_chi0_wk(g0_wk_big, nw=int(factor*nw + 1))

print g0_wk.data.shape
print chi00_wk.data.shape

U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(1, U, 0.0, 0.0, 0.0)

chi_s = solve_rpa_PH(chi00_wk, +U_s)
chi_c = solve_rpa_PH(chi00_wk, -U_c) # Minus for correct charge rpa equation

gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)

gamma_dyn, gamma_const = split_into_dynamic_and_constant(gamma)

# -- Test the Eliashberg equation with the new gamma

delta = g0_wk.copy()
delta.data[:] = 1.0

F_wk = eliashberg_g_delta_g_product(g0_wk, delta)
next_delta = eliashberg_product(gamma, g0_wk, delta)
next_delta_fft = eliashberg_product_fft(gamma, g0_wk, delta)


k_point = Idx(0, 0, 0)

print(np.allclose(next_delta[:, k_point].data, next_delta_fft[:, k_point].data))
print(np.max(np.abs((next_delta[:, k_point].data - next_delta_fft[:, k_point].data))))


from pytriqs.plot.mpl_interface import oplot, plt
subp = [3, 3, 1]
plt.figure(figsize=(9, 6))

plt.subplot(*subp); subp[-1] += 1
oplot(g0_wk[:, k_point], label='d')
plt.title('g0_wk')

plt.subplot(*subp); subp[-1] += 1
oplot(chi00_wk[:, k_point], label='d')
plt.title('chi00_wk')

plt.subplot(*subp); subp[-1] += 1
oplot(gamma[:, k_point], label='d')
plt.title('gamma')

plt.subplot(*subp); subp[-1] += 1
oplot(gamma_dyn[:, k_point], label='d')
plt.title('gamma_dyn')

plt.subplot(*subp); subp[-1] += 1
oplot(F_wk[:, k_point])
plt.title('F_wk')

plt.subplot(*subp); subp[-1] += 1
oplot(next_delta[:, k_point], label='d')
plt.title('next_delta')

plt.subplot(*subp); subp[-1] += 1
oplot(next_delta_fft[:, k_point], label='dfft')
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
