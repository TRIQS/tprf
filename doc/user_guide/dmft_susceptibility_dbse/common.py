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

import copy
import itertools

import numpy as np

from datetime import datetime

from h5 import HDFArchive
import triqs.utility.mpi as mpi

from triqs.gf import Gf, MeshImFreq, Fourier, BlockGf, inverse

from triqs.operators import c as c_operator
from triqs.operators import Operator, dagger
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.operators.util.U_matrix import U_matrix_kanamori
from triqs.operators.util.op_struct import set_operator_structure

import triqs_cthyb

from triqs_tprf.lattice import lattice_dyson_g_w
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections
from triqs_tprf.utilities import BlockGf_data

from triqs.gf.meshes import MeshDLR
from triqs_tprf.fitdlr import fitdlr
from triqs_tprf.fitdlr import BlockSymmetrizer


def setup_dmft_calculation(p):

    p = copy.deepcopy(p)
    p.iter = 0

    # -- Block structure of GF
    
    p.num_orbitals = 3
    p.spin_names = ('up','do')
    p.orb_names = list(range(p.num_orbitals))

    p.init.gf_struct = set_operator_structure(p.spin_names, p.num_orbitals, off_diag=False)
    mpi.report(f'p.init.gf_struct = \n{p.init.gf_struct}')

    p.fundamental_operators = [
        c_operator(f'{s}_{o}', 0) for s, o in
        itertools.product(p.spin_names, p.orb_names)]
    
    # -- Local Hamiltonian

    p.J = p.U * p.J_over_U
    KanMat1, KanMat2 = U_matrix_kanamori(p.num_orbitals, p.U, p.J)
    p.solve.h_int = h_int_kanamori(p.spin_names, p.num_orbitals, KanMat1, KanMat2, p.J, False)

    # -- Sr2RuO4 Wannier90 model from GPAW
    
    from tight_binding_model import tight_binding_model
    H = tight_binding_model()

    p.kmesh = H.get_kmesh(n_k = (p.n_k, p.n_k, p.n_k))
    p.e_k = H.fourier(p.kmesh)
    
    # -- Initial zero guess for the self-energy
    
    p.sigma_w = Gf(mesh=MeshImFreq(p.init.beta, 'Fermion', p.init.n_iw), target_shape=[6, 6])
    p.sigma_w.zero()
    p.sigma_w << +p.mu * np.eye(2*p.num_orbitals)
    
    return p


def solve_self_consistent_dmft_fix_N(p, verbose=False, filename=None):

    p0 = copy.deepcopy(p)
    ps = []

    def target_function(mu, p0, ps):
        mpi.report(f'--> solve_self_consistent_dmft_fix_N: target_function with mu = {mu}')

        if len(ps) != 0: p0 = copy.deepcopy(ps[-1])
        p0.mu = mu
        
        ps += solve_self_consistent_dmft(p0, verbose=False)

        if filename:
            if mpi.is_master_node():
                timestamp = f'{datetime.now().timestamp()}'
                tmp_filename = filename.replace('data_', 'data_tmp_').replace('.h5', '_time_'+timestamp+'.h5')
                print(f'--> Writing: {tmp_filename}')
                with HDFArchive(tmp_filename, 'w') as a:
                    #a['ps'] = ParameterCollections(ps)
                    a['ps'] = ps
            
        p = ps[-1]
        mpi.report(f'--> solve_self_consistent_dmft_fix_N: target_function ' + \
                   f'iter = {p.iter} mu = {p.mu}')
        mpi.report(f'--> solve_self_consistent_dmft_fix_N: target_function ' + \
                   f'N = {p.N}, N_target = {p.N_target}, N - N_target = {p.N - p.N_target}')
        return p.N - p.N_target

    from scipy.optimize import root_scalar

    sol = root_scalar(
        target_function, bracket=[p.mu_min, p.mu_max], args=(p0, ps),
        maxiter=p.mu_iter_max, xtol=p.N_tol)

    mpi.report('='*72)
    mpi.report(f'mu = {sol.root}')
    mpi.report('='*72)

    mu = sol.root
    p0 = copy.deepcopy(ps[-1])
    p0.mu = mu
    ps += solve_self_consistent_dmft(p0, verbose=False)

    return ps


def solve_self_consistent_dmft(p, verbose=False):

    ps = []
    for dmft_iter in range(p.sc_iter_max):
        mpi.report(f'--> DMFT Iteration: {p.iter}')
        p = dmft_self_consistent_step(p, verbose=verbose)
        ps.append(p)
        mpi.report(f'--> DMFT Convergence: dG = {p.dG:2.2E}')
        mpi.report(f'--> Density: N = {p.N:2.2E}')
        mpi.report(f'--> Chempot: mu = {p.mu:2.2E}')
        mpi.report(f'--> Densmat: rho_up = \n{p.rho.real[:3,:3]}')
        mpi.report(f'--> Densmat: rho_do = \n{p.rho.real[3:,3:]}')
        if p.dG < p.G_tol: break

    if p.dG > p.G_tol: mpi.report('--> Warning: DMFT Not converged!')
    else: mpi.report(f'--> DMFT Converged: dG = {p.dG:2.2E}')
    return ps


def dmft_self_consistent_step(p, verbose=False):

    p = copy.deepcopy(p)
    p.iter += 1

    p.g_w = lattice_dyson_g_w(p.mu, p.e_k, p.sigma_w)
    
    p.rho = p.g_w.density()
    p.N = np.sum(np.diagonal(p.rho)).real

    # -- Solve Impurity problem

    p.g0_w = p.g_w.copy()
    p.g0_w << inverse(inverse(p.g_w) + p.sigma_w)
    
    cthyb = triqs_cthyb.Solver(**p.init.dict())

    BlockGf_from_Gf_matrix_valued(cthyb.G0_iw, p.g0_w)

    solve = copy.deepcopy(p.solve)
    solve.n_cycles = int(np.round(solve.n_cycles / mpi.size))

    cthyb.solve(**solve.dict())

    p.G0_w = cthyb.G0_iw.copy()
    p.G_tau_raw = cthyb.G_tau
    p.Delta_infty = cthyb.Delta_infty
    
    G_tau = fit_dlr(p.G_tau_raw, p)

    p.dG = np.max(np.abs(BlockGf_data(p.G_tau - G_tau))) \
        if hasattr(p, 'G_tau') else float('nan')

    p.G_tau = G_tau

    p.G_w = p.G0_w.copy()
    p.G_w << Fourier(p.G_tau)

    p.Sigma_w = p.G0_w.copy()    
    p.Sigma_w << inverse(p.G0_w) - inverse(p.G_w)
    
    Gf_matrix_valued_from_BlockGf(p.sigma_w, p.Sigma_w)
    
    return p


def fit_dlr(G_tau, p):

    mpi.report(f'--> Fit DLR')

    G_tau = G_tau.copy()
    cmesh = MeshDLR(p.init.beta, 'Fermion', w_max=p.w_max, eps=p.eps)
    print(cmesh)

    block_mat = np.diag([1, 1, 2, 1, 1, 2])
    sym = BlockSymmetrizer(len(cmesh), block_mat)

    opt = dict(
        discontinuity=True,
        density=True,
        realvalued=True,
        ftol=1e-6,
        symmetrizer=sym,
        verbose=mpi.is_master_node(),
        )

    G_tau_mat = Gf(mesh=G_tau.mesh, target_shape=[6, 6])
    Gf_matrix_valued_from_BlockGf(G_tau_mat, G_tau)

    H_delta = Operator()
    for i in range(6):
        op = p.fundamental_operators[i]
        H_delta += p.Delta_infty[i][0,0].real * dagger(op) * op 

    H = p.solve.h_int + H_delta

    G_c_mat, sol = fitdlr(cmesh, G_tau_mat, H, p.fundamental_operators, **opt)

    from triqs.gf.gf_factories import make_gf_imtime
    G_mat_fit = make_gf_imtime(G_c_mat, len(G_tau_mat.mesh))

    BlockGf_from_Gf_matrix_valued(G_tau, G_mat_fit)
    return G_tau


def get_matrix_block_pairs():
    blocks = [ f'{s}_{o}' for s, o in itertools.product( ['up', 'do'], range(3) ) ]
    matrix_block_pairs = [ (idx, bidx) for idx, bidx in enumerate(blocks) ]
    return matrix_block_pairs


def Gf_matrix_valued_from_BlockGf(G_mat, G_block):
    for idx, bidx in get_matrix_block_pairs():
        G_mat[idx, idx] << G_block[bidx][0, 0]
    return G_mat


def BlockGf_from_Gf_matrix_valued(G_block, G_mat):
    for idx, bidx in get_matrix_block_pairs():
        G_block[bidx][0,0] << G_mat[idx, idx]
    return G_block

    
