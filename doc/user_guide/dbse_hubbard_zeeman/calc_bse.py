################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by Hugo U. R. Strand
# Author: Hugo U. R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import time
import numpy as np
import itertools

from h5 import HDFArchive

from triqs.gfs import Gf, Fourier
from triqs.gfs import make_gf_from_fourier

from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.rpa_tensor import kanamori_quartic_tensor

from triqs_tprf.bse import impurity_irreducible_vertex_Gamma
from triqs_tprf.dbse import solve_lattice_dbse
from triqs_tprf.dbse import impurity_reducible_vertex_F
from triqs_tprf.bse import solve_lattice_bse

from triqs_tprf.utilities import G2_loc_fixed_fermionic_window_python

from w2dyn_cthyb.converters_worm import p2_from_w2dyn_P2_worm_components
from w2dyn_cthyb.converters_worm import p2_remove_disconnected
from w2dyn_cthyb.converters_worm import p3_from_w2dyn_P3_worm_components
from w2dyn_cthyb.converters_worm import p3_w2dyn_to_triqs_freq_shift_alt
from w2dyn_cthyb.converters_worm import L_from_g3
from w2dyn_cthyb.converters_worm import g2_from_w2dyn_G2_worm_components

from triqs.gfs import Gf, MeshProduct, Idx, MeshImFreq
from triqs.gfs import inverse
from triqs_tprf.lattice import fourier_wk_to_wr, chi0r_from_gr_PH, chi0q_from_chi0r, chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH

from triqs_tprf.linalg import product_PH, inverse_PH, identity_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH

from triqs_tprf.bse import get_chi0_nk_at_specific_w
from triqs_tprf.lattice_utils import add_fake_bosonic_mesh

import numpy as np

def load_h5(filename):
    print(f'--> Loading: {filename}')
    with HDFArchive(filename, 'r') as a:
        p = a['p']
    return p

filename_sc  = 'data_sc.h5'
filename_chi = 'data_susc.h5'
filename_tri = 'data_triangle.h5'
filename_g2  = 'data_g2.h5'

print(f'--> Loading: {filename_sc}')
with HDFArchive(filename_sc, 'r') as a:
    p = a['ps'][-1]

for i, j in itertools.product(range(2), repeat=2):
    if i != j:
        p.g_w[i, j] = 0.


p.num_orbitals=1


# Impurity susceptibility (one frequency)
p_chi = load_h5(filename_chi)
p2 = p2_from_w2dyn_P2_worm_components(p_chi.GF_worm_components, p.num_orbitals)
p.g_tau = make_gf_from_fourier(p.g_w)
p.chi_imp_w = p2_remove_disconnected(p2, p.g_tau)

# "Triangle" impurity two-particle Green's function (two frequencies)
p_tri = load_h5(filename_tri)
p3 = p3_from_w2dyn_P3_worm_components(p_tri.GF_worm_components, p.num_orbitals)
p3 = p3_w2dyn_to_triqs_freq_shift_alt(p3)
p.L_wn = L_from_g3(p3, p.g_w) # remove disconnected and amputate

# "Square" impurity two-particle Green's function (three frequencies)
p_g2 = load_h5(filename_g2)
p.g2_wnn = g2_from_w2dyn_G2_worm_components(
    p_g2.G2_worm_components, p.num_orbitals)

# Lattice dispersion and Green's function
B_field = np.diag([+p.B, -p.B,]) if hasattr(p, 'B') else zero
g_wk = lattice_dyson_g_wk(mu=p.mu, e_k=p.e_k, sigma_w=p.sigma_w- B_field)

for nwf in [30,28,26,24,22,20,18]:
    print('='*72)
    print(f'nwf = {nwf}', flush=True)
    p.nwf = nwf
    g2_wnn = G2_loc_fixed_fermionic_window_python(p.g2_wnn, nwf=p.nwf)

    # G2 -> F
    p.F_wnn = impurity_reducible_vertex_F(p.g_w, g2_wnn)
    p.chi_kw_dbse = solve_lattice_dbse(g_wk, p.F_wnn, p.L_wn, p.chi_imp_w)

    print('--> BSE (for reference)')
    Gamma_wnn = impurity_irreducible_vertex_Gamma(p.g_w, g2_wnn)
    p.chi_kw_bse, p.chi0_kw = solve_lattice_bse(g_wk, Gamma_wnn)

    filename_out = f'data_bse_nwf_{nwf:03d}.h5'
    print(f'--> Saving: {filename_out}')
    with HDFArchive(filename_out, 'w') as a:
        a['p'] = p
