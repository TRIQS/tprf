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

from triqs.gf import Gf, Fourier
from triqs.gf import make_gf_from_fourier

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

# DBSEP

# Implementation of dual Bethe-Salpeter equation for the polarization (DBSEP) according to [1] Krien, PHYSICAL REVIEW B 99, 235106 (2019)
# For the TPRF implementation of the dual Bethe-Salpeter equation for the susceptibility, see [2] van Loon and Strand, PHYSICAL REVIEW B 109, 155157 (2024)
# Erik van Loon, Lund University, 2025
# Experimental! Make sure to verify your results against regular BSE or DBSE

# A summary of the method is as follows:
# The dual Bethe-Salpeter equation (DBSE) is already more efficient than the BSE, since it uses Gdual=G-gloc instead of G,
# which decays as 1/nu^2 instead of 1/nu and therefore leads to a quickly decaying bubble and 1/nwf convergence of the DBSE,
# compared to 1/nwf for the ordinary BSE.
# The DBSE corresponds to a resummation of processes with local fermionic propagators, and this resummation leads to the desired replacement G->Gdual
# In the same spirit, a remaining bottleneck for the convergence is the fact that the vertex F in the DBSE is asymptotically constant.
# Krien's formula (here called DBSEP) is a resummation of the susceptibility in terms of processes with reducible local interactions. 
# It replaces F by Firr, where Firr has the property that it asymptotically decays instead of being constant. 
# Thus, all vertex corrections in the DBSEP converge with a higher power of 1/nwf than those in the DBSE
#
# Note: in the DBSEP, processes which don't involve Firr, only Lirr, also contribute asympotically, but they are computationally simpler. 
# The DBSEP has better formal convergence only if those processes are taken into account using a sufficiently larger frequency box, which is not implemented here.
# Instead, the current implementation is expected to have the same power of nwf but better prefactor.

from triqs.gf import Gf, MeshProduct, Idx, MeshImFreq
from triqs.gf import inverse
from triqs_tprf.lattice import fourier_wk_to_wr, chi0r_from_gr_PH, chi0q_from_chi0r, chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH

from triqs_tprf.linalg import product_PH, inverse_PH, identity_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH

from triqs_tprf.bse import get_chi0_nk_at_specific_w
from triqs_tprf.lattice_utils import add_fake_bosonic_mesh

import numpy as np

from linalg_utils import ChannelOrder
from linalg_utils import matrix_from_tensor, tensor_from_matrix
from linalg_utils import gf_matrix_from_tensor, gf_tensor_from_matrix

def impurity_polarization(chi_w,U_tensor):

    chi_w_mat = gf_matrix_from_tensor(chi_w)
    U_mat = matrix_from_tensor(U_tensor)
    
    pi_w_mat = chi_w_mat.copy()
    W_w_mat = chi_w_mat.copy()

    # Definition above Eq 8 in [1]:
    # chi = pi / (1 - U pi), so chi - chi U pi = pi. Due to geometric series, also chi - pi U chi = pi, sochi = pi (1+ U chi) and finally:
    # pi = chi / (1+U chi)
    # W  = U / (1 - pi U )

    pi_w_mat << chi_w_mat * inverse(1. + U_mat * chi_w_mat )
    W_w_mat << U_mat * inverse(1. - pi_w_mat*U_mat)

    pi_w = gf_tensor_from_matrix(pi_w_mat)
    W_w = gf_tensor_from_matrix(W_w_mat)

    return pi_w, W_w

def pi_kw_to_chi_kw(pi_kw, U_tensor):
    
    U_mat = matrix_from_tensor(U_tensor)
    pi_kw_mat = gf_matrix_from_tensor(pi_kw)
    chi_kw_mat= pi_kw_mat.copy()

    # chi = pi / (1 - U pi)
    for k in pi_kw.mesh[0]:
        chi_kw_mat[k,:] << pi_kw_mat[k,:] * inverse(1. - U_mat*pi_kw_mat[k,:])

    return gf_tensor_from_matrix(chi_kw_mat)

def irreducible_L(L_wn, U_tensor, pi_w):
    U_mat = matrix_from_tensor(U_tensor)
    pi_w_mat = gf_matrix_from_tensor(pi_w)
    L_wn_mat = gf_matrix_from_tensor(L_wn)
    Lirr_wn_mat = L_wn_mat.copy()

    # L = Lirr / (1-U pi), according to Eq 8 of [1]. Solving for Lirr gives
    # L - L U pi = Lirr
    # Note: TPRF uses L with boson on the right, so it is L U pi and not pi U L
    for n in L_wn.mesh[1]:
        Lirr_wn_mat[:,n] = L_wn_mat[:,n] - L_wn_mat[:,n] * U_mat * pi_w_mat

    return gf_tensor_from_matrix(Lirr_wn_mat)

def irreducible_F(F_wnn,W_w, L_wn):
    # Make irreducible vertex Firr from reducible version F
    # Eq 9 of [1]

    Firr_wnn = F_wnn.copy()
    for w,n1,n2 in F_wnn.mesh:
        Firr_wnn[w,n1,n2] = F_wnn[w,n1,n2] - np.einsum('abfe,efgh,dcgh->abcd', L_wn[Idx(w.index),Idx(n1.index)], W_w[Idx(w.index)], L_wn[Idx(-w.index),Idx(-n2.index-1) ].conj(), optimize=True ) 
        # Mirroring implemented using symmetry of L
        # Note the conventions of the L vertex in Refs [1,2]
        # TPRF's DBSE uses L as a vertex with the bosonic line on the right (Fig 2 of [2])
        # Eq 9 / Fig 1b of [1] have two L vertex, with the bosonic line on the right and left, respectively.
        # Eq 16 of [1] can be used to to map between right-sided and left-sided L-vertices, using complex conjugation
    return Firr_wnn


###

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

# Remove small (1e-6) off diagonal terms in e_k and g_w by hand
#e_loc = np.sum(p.e_k.data, axis=0).real / p.e_k.data.shape[0]
#e_loc -= np.diag(np.diag(e_loc))
#p.e_k.data[:] -= e_loc[None, ...]

for i, j in itertools.product(range(2), repeat=2):
    if i != j:
        p.g_w[i, j] = 0.


# Interaction 
# Note: this TPRF implementation does the Fierz ambiguity nicely, so that Firr=0 asymptotically
p.U_mat = kanamori_quartic_tensor(1, p.U, p.U, 0, 0)
p.num_orbitals=1


# Impurity susceptibility (one frequency)
p_chi = load_h5(filename_chi)
p2 = p2_from_w2dyn_P2_worm_components(p_chi.GF_worm_components, p.num_orbitals)
p.g_tau = make_gf_from_fourier(p.g_w)
p.chi_imp_w = p2_remove_disconnected(p2, p.g_tau)
# Impurity polarization/screened interaction
p.pi_imp_w, p.W_w = impurity_polarization( p.chi_imp_w, p.U_mat)

# "Triangle" impurity two-particle Green's function (two frequencies)
p_tri = load_h5(filename_tri)
p3 = p3_from_w2dyn_P3_worm_components(p_tri.GF_worm_components, p.num_orbitals)
p3 = p3_w2dyn_to_triqs_freq_shift_alt(p3)
p.L_wn = L_from_g3(p3, p.g_w) # remove disconnected and amputate
# U-irreducible part
p.Lirr_wn = irreducible_L(p.L_wn, p.U_mat, p.pi_imp_w)

# "Square" impurity two-particle Green's function (three frequencies)
p_g2 = load_h5(filename_g2)
p.g2_wnn = g2_from_w2dyn_G2_worm_components(
    p_g2.G2_worm_components, p.num_orbitals)

# Lattice dispersion and Green's function
B_field = np.diag([+p.B, -p.B,]) if hasattr(p, 'B') else zero
g_wk = lattice_dyson_g_wk(mu=p.mu, e_k=p.e_k, sigma_w=p.sigma_w- B_field)

# DBSEP calculations for varying frequency window

for nwf in [30,28,26,24,22,20,18]:
#for nwf in [4,]:
    print('='*72)
    print(f'nwf = {nwf}', flush=True)
    p.nwf = nwf
    g2_wnn = G2_loc_fixed_fermionic_window_python(p.g2_wnn, nwf=p.nwf)

    # G2 -> F
    p.F_wnn = impurity_reducible_vertex_F(p.g_w, g2_wnn)
    # U-irreducible part
    p.Firr_wnn = irreducible_F(p.F_wnn, p.W_w, p.Lirr_wn)

    p.chi_kw_dbse = solve_lattice_dbse(g_wk, p.F_wnn, p.L_wn, p.chi_imp_w)
    p.pi_kw = solve_lattice_dbse(g_wk, p.Firr_wnn, p.Lirr_wn, p.pi_imp_w) # Same function call, different arguments
    p.chi_kw_dbsep = pi_kw_to_chi_kw(p.pi_kw, p.U_mat)

    print('--> BSE (for reference)')
    chi_wnn = chi_from_gg2_PH(p.g_w, g2_wnn)
    Gamma_wnn = impurity_irreducible_vertex_Gamma(p.g_w, chi_wnn)
    p.chi_kw_bse, p.chi0_kw = solve_lattice_bse(g_wk, Gamma_wnn)

    filename_out = f'data_dbsep_nwf_{nwf:03d}.h5'
    print(f'--> Saving: {filename_out}')
    with HDFArchive(filename_out, 'w') as a:
        a['p'] = p
   
