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
import numpy as np

import triqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs.gf import Gf, MeshImFreq, Fourier, LegendreToMatsubara, BlockGf, inverse, Idx

import triqs_cthyb

from triqs_tprf.lattice import lattice_dyson_g_w
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections
from triqs_tprf.utilities import BlockGf_data

def setup_dmft_calculation(p):

    p = copy.deepcopy(p)
    p.iter = 0

    # -- Local Hubbard interaction
    from triqs.operators import n
    p.solve.h_int = p.U*n('up', 0)*n('do', 0) - 0.5*p.U*(n('up', 0) + n('do', 0))

    # -- 2D square lattice w. nearest neighbour hopping t
    from triqs_tprf.tight_binding import TBLattice
    T = -p.t * np.eye(2)
    H = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        orbital_positions = [(0,0,0)]*2,
        orbital_names = ['up', 'do'],
        hopping = {(0, +1) : T, (0, -1) : T, (+1, 0) : T, (-1, 0) : T})
    
    kmesh = H.get_kmesh(n_k = (p.n_k, p.n_k, 1))
    p.e_k = H.fourier(kmesh)

    # -- Initial zero guess for the self-energy    
    p.sigma_w = Gf(mesh=MeshImFreq(p.init.beta, 'Fermion', p.init.n_iw), target_shape=[2, 2])
    p.sigma_w.zero()

    return p

def solve_self_consistent_dmft(p):

    ps = []
    for dmft_iter in range(p.n_iter):
        mpi.report('--> DMFT Iteration: {:d}'.format(p.iter))
        p = dmft_self_consistent_step(p)
        ps.append(p)
        mpi.report('--> DMFT Convergence: dG_l = {:f}'.format(p.dG_l))
        if p.dG_l < p.G_l_tol: break

    if dmft_iter >= p.n_iter - 1: mpi.report('--> Warning: DMFT Not converged!')
    else: mpi.report('--> DMFT Converged: dG_l = {:f}'.format(p.dG_l))
    return ps

def dmft_self_consistent_step(p):

    p = copy.deepcopy(p)
    p.iter += 1

    p.g_w = lattice_dyson_g_w(p.mu, p.e_k, p.sigma_w - np.diag([p.B, -p.B]))
    p.g0_w = p.g_w.copy()
    p.g0_w << inverse(inverse(p.g_w) + p.sigma_w)

    cthyb = triqs_cthyb.Solver(**p.init.dict())

    # -- set impurity from lattice
    cthyb.G0_iw['up'][0, 0] << p.g0_w[0, 0]
    cthyb.G0_iw['do'][0, 0] << p.g0_w[1, 1]

    cthyb.solve(**p.solve.dict())

    p.dG_l = np.max(np.abs(BlockGf_data(cthyb.G_l-p.G_l))) if hasattr(p,'G_l') else float('nan')
    p.G_l = cthyb.G_l

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
    M_old = p.M if hasattr(p, 'M') else float('nan')
    p.M = 0.5*(p.rho[0, 0] - p.rho[1, 1])
    p.dM = np.abs(p.M - M_old)
    
    return p
