################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
# Author: H. U.R. Strand
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

import copy
import glob
import itertools
import numpy as np

from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs.gf import Gf, MeshImFreq, Fourier, LegendreToMatsubara, BlockGf, inverse, Idx
from triqs.gf import MeshDLRImFreq, make_gf_dlr, fit_gf_dlr, make_gf_dlr_imfreq, make_gf_dlr_imtime
from triqs.gf import make_gf_imfreq, make_gf_imtime

from triqs.operators import c as c_operator
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.operators.util.U_matrix import U_matrix_kanamori

import triqs_cthyb

from triqs_tprf.lattice import dlr_on_imfreq
from triqs_tprf.lattice import lattice_dyson_g_w
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections
from triqs_tprf.utilities import BlockGf_data


def setup_dmft_calculation(p):

    p = copy.deepcopy(p)
    p.iter = 0

    p.spin_names = ('up','do')
    p.orb_names = list(range(p.num_orbitals))

    p.num_spin_orb = p.num_orbitals * 2
    p.spin_orb_names = [f'{s}_{o}' for s, o in
        itertools.product(p.spin_names, p.orb_names)]
    
    p.fundamental_operators = [
        c_operator(name, 0) for name in p.spin_orb_names]

    p.init.gf_struct = [ (name, 1) for name in p.spin_orb_names ]
    
    # -- Local Hubbard interaction
    if p.num_orbitals == 1:
        from triqs.operators import n
        p.solve.h_int = p.U*n('up_0', 0)*n('do_0', 0)
    else:
        p.J = p.U * p.J_over_U
        KanMat1, KanMat2 = U_matrix_kanamori(p.num_orbitals, p.U, p.J)
        p.solve.h_int = h_int_kanamori(p.spin_names, p.num_orbitals, KanMat1, KanMat2, p.J, False)

    # -- 2D square lattice w. nearest neighbour hopping t
    from triqs_tprf.tight_binding import TBLattice
    T0 = np.kron(np.diag([p.B, -p.B]), np.eye(p.num_orbitals))
    T = -p.t * np.eye(p.num_spin_orb)
    H = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        orbital_positions = [(0,0,0)]*p.num_spin_orb,
        orbital_names = p.spin_orb_names,
        hopping = {(0, 0) : T0,
                   (0, +1) : T, (0, -1) : T, (+1, 0) : T, (-1, 0) : T})
    
    kmesh = H.get_kmesh(n_k = (p.n_k, p.n_k, 1))
    p.e_k = H.fourier(kmesh)

    # -- Initial zero guess for the self-energy
    p.wmesh = MeshDLRImFreq(p.init.beta, 'Fermion', p.w_max, p.eps)

    p.sigma_w = Gf(mesh=p.wmesh, target_shape=[p.num_spin_orb, p.num_spin_orb])
    p.sigma_w.zero()
    p.sigma_w.data[:] = +p.U/2 * np.eye(p.num_spin_orb)[None, ...]

    sigma_c = make_gf_dlr(p.sigma_w)
    sigma_t = make_gf_dlr_imtime(sigma_c)
    
    p.cmesh = sigma_c.mesh
    p.tmesh = sigma_t.mesh

    return p

def solve_self_consistent_dmft(p):

    ps = []
    for dmft_iter in range(p.n_iter):
        mpi.report('--> DMFT Iteration: {:d}'.format(p.iter))
        p = dmft_self_consistent_step(p)
        ps.append(p)
        mpi.report('--> DMFT Convergence: dG = {:2.2E}'.format(p.dG))
        if p.dG < p.G_tol: break

    if dmft_iter >= p.n_iter - 1: mpi.report('--> Warning: DMFT Not converged!')
    else: mpi.report('--> DMFT Converged: dG = {:2.2E}'.format(p.dG))
    return ps

def dmft_self_consistent_step(p):

    p = copy.deepcopy(p)
    p.iter += 1

    p.g_w = lattice_dyson_g_w(p.mu, p.e_k, p.sigma_w)

    p.g_c = make_gf_dlr(p.g_w)
    p.g0_w = p.g_w.copy()
    p.g0_w << inverse(inverse(p.g_w) + p.sigma_w)
    p.g0_c = make_gf_dlr(p.g0_w)

    cthyb = triqs_cthyb.Solver(**p.init.dict())

    # -- set impurity from lattice
    wmesh = cthyb.G0_iw.mesh
    g0_w = dlr_on_imfreq(p.g0_c, wmesh)
    BlockGf_from_Gf_matrix_valued(cthyb.G0_iw, g0_w)

    cthyb.solve(**p.solve.dict())

    p.G_tau_raw = cthyb.G_tau.copy()
    p.G0_w, p.G_w, p.Sigma_w = cthyb.G0_iw.copy(), cthyb.G_iw.copy(), cthyb.Sigma_iw.copy()

    p.orbital_occupations = cthyb.orbital_occupations
    p.G_moments = cthyb.G_moments
    p.Sigma_moments = cthyb.Sigma_moments
    p.Delta_infty = cthyb.Delta_infty

    p.rho = symmetrize_mat(block_dict_to_mat(p.orbital_occupations, p.spin_orb_names))
    p.G2 = symmetrize_mat(moment_block_dict_to_mat(p.G_moments, 2, p.spin_orb_names))

    p.S0 = symmetrize_mat(moment_block_dict_to_mat(p.Sigma_moments, 0, p.spin_orb_names))
    p.S1 = symmetrize_mat(moment_block_dict_to_mat(p.Sigma_moments, 1, p.spin_orb_names))
    p.D0 = symmetrize_mat(np.diag(np.array(p.Delta_infty).flatten()))

    p.G3 = p.S1 + (p.D0 + p.S0) @ (p.D0 + p.S0)

    # -- Fit

    if False:
        p.G_c = fit_gf_dlr(p.G_tau_raw, p.w_max, p.eps)
        p.G_tau = make_gf_imtime(p.G_c, n_tau=p.init.n_tau)
        p.G_w = make_gf_dlr_imfreq(p.G_c)

        p.g_w_cthyb = p.g_w.copy()
        Gf_matrix_valued_from_BlockGf(p.g_w_cthyb, p.G_w)
        
    else:
        p.G_tau, p.g_c = fit_dlr(p.G_tau_raw, p)
        p.g_w_cthyb = make_gf_dlr_imfreq(p.g_c)
        p.g_w_linear = make_gf_imfreq(p.g_c, n_iw=len(p.G_w.mesh))

    # -- Difference between iterations
    
    p.dG = np.max(np.abs(p.g_w.data - p.g_w_cthyb.data))
    
    # -- set lattice from impurity
    
    p.sigma_w << inverse(p.g0_w) - inverse(p.g_w_cthyb)

    tmp = p.sigma_w.copy()
    tmp.data[:] -= p.S0[None, ...]
    p.sigma_c = make_gf_dlr(tmp)
    p.sigma_w_linear = make_gf_imfreq(p.sigma_c, n_iw=len(p.G_w.mesh)) + p.S0
    
    # -- local observables
    p.rho = p.g_c.density()
    M_old = p.M if hasattr(p, 'M') else float('nan')
    p.M = 0.5*(p.rho[0, 0] - p.rho[1, 1])
    p.dM = np.abs(p.M - M_old)
    
    return p


def fit_dlr(G_tau, p):

    #from triqs_tprf.fitdlr import fitdlr
    #from triqs_tprf.fitdlr import BlockSymmetrizer

    from fitdlr import fitdlr
    from fitdlr import BlockSymmetrizer

    mpi.report(f'--> Fit DLR')

    G_tau = G_tau.copy()
    cmesh = p.cmesh

    block_mat = np.diag([1, 2])
    sym = BlockSymmetrizer(len(cmesh), block_mat)

    opt = dict(
        discontinuity=True,
        #density=True,
        density=False,
        realvalued=True,
        positivity=True,
        ftol=1e-7,
        symmetrizer=sym,
        verbose=mpi.is_master_node(),
        rho=p.rho,
        G2=p.G2,
        #G3=p.G3,
        )

    G_tau_mat = Gf(mesh=G_tau.mesh, target_shape=[p.num_spin_orb, p.num_spin_orb])
    Gf_matrix_valued_from_BlockGf(G_tau_mat, G_tau)

    from triqs.operators import Operator, dagger
    H_delta = Operator()
    for i in range(len(p.fundamental_operators)):
        op = p.fundamental_operators[i]
        H_delta += p.Delta_infty[i][0,0].real * dagger(op) * op 

    H = p.solve.h_int + H_delta

    G_c_mat, sol = fitdlr(cmesh, G_tau_mat, H, p.fundamental_operators, **opt)

    from triqs.gf.gf_factories import make_gf_imtime
    G_mat_fit = make_gf_imtime(G_c_mat, len(G_tau_mat.mesh))

    BlockGf_from_Gf_matrix_valued(G_tau, G_mat_fit)
    return G_tau, G_c_mat


def get_matrix_block_pairs():
    blocks = [ f'{s}_{o}' for s, o in itertools.product( ['up', 'do'], range(num_orbitals) ) ]
    matrix_block_pairs = [ (idx, bidx) for idx, bidx in enumerate(blocks) ]
    return matrix_block_pairs


def Gf_matrix_valued_from_BlockGf(G_mat, G_block):
    for idx, bidx in enumerate(G_block.indices):
        G_mat[idx, idx] << G_block[bidx][0, 0]
    return G_mat


def BlockGf_from_Gf_matrix_valued(G_block, G_mat):
    for idx, bidx in enumerate(G_block.indices):
        G_block[bidx][0,0] << G_mat[idx, idx]
    return G_block


def block_dict_to_mat(bd, spin_orb_names):
    nso = len(spin_orb_names)
    mat = np.zeros([nso, nso], dtype=complex)
    for idx, bidx in enumerate(spin_orb_names):
        mat[idx, idx] = bd[bidx]
    return mat


def moment_block_dict_to_mat(bd, moment, spin_orb_names):
    nso = len(spin_orb_names)
    mat = np.zeros([nso, nso], dtype=complex)
    for idx, bidx in enumerate(spin_orb_names):
        mat[idx, idx] = bd[bidx][moment]
    return mat


def symmetrize_mat(mat):
    #avg = np.mean(np.diag(mat))
    #mat = np.diag([avg]*2)
    return mat
