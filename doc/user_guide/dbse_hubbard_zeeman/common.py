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
import glob
import itertools
from datetime import datetime

import numpy as np

import triqs.utility.mpi as mpi
from h5 import HDFArchive

from triqs.gf import Gf, MeshImFreq, Fourier, LegendreToMatsubara, BlockGf, inverse, Idx, MeshProduct
from triqs.gf.tools import fit_legendre

from triqs.operators import n, c, c_dag, Operator
from triqs.operators.util.hamiltonians import h_int_kanamori, h_int_density
from triqs.operators.util.op_struct import set_operator_structure

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import lattice_dyson_g_w
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections
from triqs_tprf.utilities import BlockGf_data

from triqs_tprf.linalg import product_PH

from triqs.operators import c as c_operator
from triqs.operators.util.U_matrix import U_matrix_kanamori
from triqs.operators.util.op_struct import set_operator_structure

def setup_dmft_calculation(p):

    p = copy.deepcopy(p)
    p.iter = 0

    # Block structure of GF
    num_orbitals = 1
    spin_names = ('up','do')
    orb_names = list(range(num_orbitals))

    p.init.gf_struct = set_operator_structure(spin_names, orb_names, off_diag=True)
    mpi.report(f'p.init.gf_struct = \n{p.init.gf_struct}')
    
    # -- Local Hamiltonian

    p.solve.move_double = True

    p.solve.h_int = p.U*n('up',0)*n('do',0)

    # -- 2D square lattice w. nearest neighbour hopping t
    # -- and next nearest neighbour hopping t_prime
    
    T = -p.t * np.eye(2)
    T_prime = -p.t_prime * np.eye(2)

    H = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        orbital_positions = [(0,0,0)]*2,
        orbital_names = ['up', 'do'],
        hopping = {
            (0, +1) : T,
            (0, -1) : T,
            (+1, 0) : T,
            (-1, 0) : T,
            # -- next nearest neighbour hopping
            (-1, -1) : T_prime,
            (-1, +1) : T_prime,
            (+1, -1) : T_prime,
            (+1, +1) : T_prime,
            },
        )

    p.kmesh = H.get_kmesh(n_k = (p.n_k, p.n_k, 1))
    p.e_k = H.fourier(p.kmesh)

    
    # -- Initial zero guess for the self-energy    
    p.sigma_w = Gf(mesh=MeshImFreq(p.init.beta, 'Fermion', p.init.n_iw), target_shape=[2,2])
    p.sigma_w.zero()

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
        mpi.report(f'--> DMFT Convergence: dG_l = {p.dG_l:2.2E}')

    if p.dG_l > p.G_l_tol: mpi.report('--> Warning: DMFT Not converged!')
    return ps

def dmft_self_consistent_step(p, verbose=False):

    p = copy.deepcopy(p)
    p.iter += 1

    zero = np.zeros((2,2))
    
    B_field = np.diag([+p.B, -p.B,]) if hasattr(p, 'B') else zero
    
    p.g_w = lattice_dyson_g_w(p.mu, p.e_k, p.sigma_w - B_field)
    p.g0_w = p.g_w.copy()
    p.g0_w << inverse(inverse(p.g_w) + p.sigma_w)

    import triqs_cthyb
    cthyb = triqs_cthyb.Solver(**p.init.dict())

    # -- set impurity from lattice
    cthyb.G0_iw['up'][0, 0] << p.g0_w[0, 0]
    cthyb.G0_iw['do'][0, 0] << p.g0_w[1, 1]

    solve = copy.deepcopy(p.solve)
    solve.n_cycles = int(np.round(solve.n_cycles / mpi.size))
    cthyb.solve(**solve.dict())

    G_l = fit_legendre(cthyb.G_tau, order=p.n_l)
    
    from triqs.gf import enforce_discontinuity
    for bidx, g_l in G_l: enforce_discontinuity(g_l, np.eye(g_l.target_shape[0]))
    
    p.dG_l = np.max(np.abs(BlockGf_data(G_l-p.G_l))) if hasattr(p,'G_l') else float('nan')
    p.G_l = G_l.copy()
        
    p.G_tau, p.G_tau_raw = cthyb.G_tau.copy(), cthyb.G_tau.copy()
    p.G0_w, p.G_w, p.Sigma_w = cthyb.G0_iw.copy(), cthyb.G0_iw.copy(), cthyb.G0_iw.copy()
    
    p.G_tau << LegendreToMatsubara(p.G_l)
    p.G_w << Fourier(p.G_tau)
    p.Sigma_w << inverse(p.G0_w) - inverse(p.G_w)

    # -- set lattice from impurity
    p.sigma_w[0, 0] << p.Sigma_w['up'][0, 0]
    p.sigma_w[1, 1] << p.Sigma_w['do'][0, 0]

    # -- local observables
    p.rho = p.g_w.density()
    p.N = np.sum(np.diagonal(p.rho)).real
    M_old = p.M if hasattr(p, 'M') else float('nan')
    p.M = 0.5*(p.rho[0, 0] - p.rho[1, 1]).real
    p.dM = np.abs(p.M - M_old)
    
    return p

def get_matrix_block_pairs():
    blocks = [ f'{s}_{o}' for s, o in itertools.product( ['up', 'do'], range(1) ) ]
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
