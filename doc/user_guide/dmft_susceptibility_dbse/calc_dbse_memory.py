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

# This script calculates the dual Bethe-Salpeter equation one bosonic frequency
# at a time, so that the memory footprint of the calculation is reduced,
# at the cost of some increase in the compute time.
# Note: this script is intended to be run with OpenMP parallelization, i.e.,
# OMP_NUM_THREADS=XXX python3 calc_dbse_memory.py
# where XXX is the number of available threads. 

import time
import numpy as np

from h5 import HDFArchive

from triqs.gf import Gf, Fourier
from triqs.gf import make_gf_from_fourier

from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.utilities import G2_loc_fixed_fermionic_window_python

from w2dyn_cthyb.converters_worm import p2_from_w2dyn_P2_worm_components
from w2dyn_cthyb.converters_worm import p2_remove_disconnected
from w2dyn_cthyb.converters_worm import p3_from_w2dyn_P3_worm_components
from w2dyn_cthyb.converters_worm import p3_w2dyn_to_triqs_freq_shift_alt
from w2dyn_cthyb.converters_worm import L_from_g3
from w2dyn_cthyb.converters_worm import g2_from_w2dyn_G2_worm_components

from tight_binding_model import tight_binding_model

def solve_lattice_dbse(g_wk, g_w, g2_wnn, L_wn, chi_imp_w):
    from triqs.gf import MeshImFreq
    from triqs.gf import Idx,MeshProduct
    from triqs_tprf.lattice import fourier_wk_to_wr
    from triqs_tprf.lattice import chi0_nr_from_gr_PH_at_specific_w
    from triqs_tprf.lattice import chi0r_from_gr_PH
    from triqs_tprf.lattice import chi0q_from_chi0r
    from triqs_tprf.lattice import chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH

    from triqs_tprf.linalg import product_PH, inverse_PH
    from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH


    bmesh = g2_wnn.mesh[0]
    fmesh = g2_wnn.mesh[1]

    assert( len(fmesh) <= len(g_wk.mesh[0]) )
    assert( len(bmesh) <= len(L_wn.mesh[0]) )
    assert( len(fmesh) <= len(L_wn.mesh[1]) )
    assert( len(bmesh) <= len(chi_imp_w.mesh) )
    
    nw = (len(bmesh) + 1) // 2
    nn = len(fmesh) // 2
    
    kmesh = g_wk.mesh[1]


    # NOTE: This can also be done for one bosonic frequency.
    F_wnn = impurity_reducible_vertex_F(p.g_w, g2_wnn)

    fake_bmesh = MeshImFreq(beta=bmesh.beta,S='Boson',n_iw=1)
    chi_kw = Gf(indices=g2_wnn.indices,mesh=MeshProduct(kmesh,bmesh))
    for w in bmesh:
        print(w,w.index)

        #Note: this can be pulled out of the loop to save cpu time at the cost of memory footprint
        # -- Remove local gf component (at r = 0)
        print('--> g_nonlocal_wr')
        g_nonlocal_wr = fourier_wk_to_wr(g_wk)
        g_nonlocal_wr[:, Idx(0, 0, 0)] = 0.

        #Loop over w to save memory
        chi0_nonlocal_nr = chi0_nr_from_gr_PH_at_specific_w(nw_index=w.index, nn=nn, g_nr=g_nonlocal_wr)
        chi0_nonlocal_wnr=Gf(indices=chi0_nonlocal_nr.indices, mesh=MeshProduct(fake_bmesh,chi0_nonlocal_nr.mesh[0],chi0_nonlocal_nr.mesh[1]))
        for fake_w,n,r in chi0_nonlocal_wnr.mesh:
            chi0_nonlocal_wnr[fake_w,n,r] = chi0_nonlocal_nr[n,r]
        del g_nonlocal_wr
        del chi0_nonlocal_nr

        print('--> chi0_nonlocal_wnk')
        chi0_nonlocal_wnk = chi0q_from_chi0r(chi0_nonlocal_wnr)

        del chi0_nonlocal_wnr

        print('--> Resize L_wn')
        L_resize_wn = Gf(indices=L_wn.indices, mesh=MeshProduct(fake_bmesh, fmesh))
        for fake_w, n in L_resize_wn.mesh:
            L_resize_wn[fake_w, n] = L_wn[Idx(w.index), Idx(n.index)]    

        F_resize_wnn = Gf(indices=F_wnn.indices, mesh=MeshProduct(fake_bmesh, fmesh, fmesh))
        for fake_w, n1, n2 in F_resize_wnn.mesh:
            F_resize_wnn[fake_w,n1,n2] = F_wnn[Idx(w.index),n1,n2]

        chi_kw_ = chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH(
        chi0_nonlocal_wnk, F_resize_wnn, L_resize_wn)

        for k,fake_w in chi_kw_.mesh:
            chi_kw[k,w] += chi_kw_[k,fake_w] + chi_imp_w[Idx(w.index)].data
        del L_resize_wn
        del F_resize_wnn
        del chi0_nonlocal_wnk
        del chi_kw_

    return chi_kw

def impurity_reducible_vertex_F(g_w, g2_wnn):
    from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH
    from triqs_tprf.linalg import product_PH, inverse_PH

    chi_wnn = chi_from_gg2_PH(g_w, g2_wnn)
    chi0_wnn = chi0_from_gg2_PH(g_w, g2_wnn)

    g2_conn_wnn = chi_wnn - chi0_wnn

    inv_chi0_wnn = inverse_PH(chi0_wnn)
    F_wnn = product_PH(inv_chi0_wnn, product_PH(g2_conn_wnn, inv_chi0_wnn))
    
    return F_wnn

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

# -- DBSE calculations for varying frequency window

for nwf in [5,]:
    print('='*72)
    print(f'nwf = {nwf}', flush=True)
    p.nwf = nwf
    g2_wnn = G2_loc_fixed_fermionic_window_python(p.g2_wnn, nwf=p.nwf)

    print('--> DBSE')
    p.chi_kw_dbse = solve_lattice_dbse(g_wk, p.g_w, g2_wnn, p.L_wn, p.chi_imp_w) 

    filename_out = f'./data/data_mem_bse_nwf_{nwf:03d}_nk_{p.n_k:03d}.h5'
    print(f'--> Saving: {filename_out}')
    with HDFArchive(filename_out, 'w') as a:
        a['p'] = p
   
