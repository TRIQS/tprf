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

from h5 import HDFArchive

from triqs.gf import Gf, Fourier
from triqs.gf import make_gf_from_fourier

from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.bse import solve_lattice_bse
from triqs_tprf.bse import impurity_irreducible_vertex_Gamma

from triqs_tprf.dbse import solve_lattice_dbse
from triqs_tprf.dbse import impurity_reducible_vertex_F

from triqs_tprf.utilities import G2_loc_fixed_fermionic_window_python

from w2dyn_cthyb.converters_worm import p2_from_w2dyn_P2_worm_components
from w2dyn_cthyb.converters_worm import p2_remove_disconnected
from w2dyn_cthyb.converters_worm import p3_from_w2dyn_P3_worm_components
from w2dyn_cthyb.converters_worm import p3_w2dyn_to_triqs_freq_shift_alt
from w2dyn_cthyb.converters_worm import L_from_g3
from w2dyn_cthyb.converters_worm import g2_from_w2dyn_G2_worm_components

from tight_binding_model import tight_binding_model


def load_h5(filename):
    print(f'--> Loading: {filename}')
    with HDFArchive(filename, 'r') as a:
        p = a['p']
    return p

filename_sc  = './data/data_sc.h5'
filename_chi = './data/data_chi.h5'
filename_tri = './data/data_tri.h5'
filename_g2  = './data/data_g2.h5'

print(f'--> Loading: {filename_sc}')
with HDFArchive(filename_sc, 'r') as a:
    p = a['ps'][-1]

# Remove small (1e-6) off diagonal terms in e_k and g_w by hand

e_loc = np.sum(p.e_k.data, axis=0).real / p.e_k.data.shape[0]
e_loc -= np.diag(np.diag(e_loc))
p.e_k.data[:] -= e_loc[None, ...]

import itertools
for i, j in itertools.product(range(6), repeat=2):
    if i != j:
        p.g_w[i, j] = 0.

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

# -- Lattice dispersion and Green's function

p.n_k = 16 # Set k-point resolution
H = tight_binding_model()
p.kmesh = H.get_kmesh(n_k = (p.n_k, p.n_k, p.n_k))
p.e_k = H.fourier(p.kmesh)
g_wk = lattice_dyson_g_wk(mu=p.mu, e_k=p.e_k, sigma_w=p.sigma_w)

# -- DBSE, BSE calculations for varying frequency window

for nwf in [40, 30, 20, 10]:
    print('='*72)
    print(f'nwf = {nwf}', flush=True)
    p.nwf = nwf
    g2_wnn = G2_loc_fixed_fermionic_window_python(p.g2_wnn, nwf=p.nwf)

    print('--> DBSE')
    p.F_wnn = impurity_reducible_vertex_F(p.g_w, g2_wnn)
    p.chi_kw_dbse = solve_lattice_dbse(g_wk, p.F_wnn, p.L_wn, p.chi_imp_w)

    print('--> BSE (for reference)')
    Gamma_wnn = impurity_irreducible_vertex_Gamma(p.g_w, g2_wnn)
    p.chi_kw_bse, p.chi0_kw = solve_lattice_bse(g_wk, Gamma_wnn)
        
    filename_out = f'./data/data_bse/data_bse_nwf_{nwf:03d}_nk_{p.n_k:03d}.h5'
    print(f'--> Saving: {filename_out}')
    with HDFArchive(filename_out, 'w') as a:
        a['p'] = p
   
