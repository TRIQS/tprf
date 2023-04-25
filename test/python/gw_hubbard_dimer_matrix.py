################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by Hugo U.R. Strand
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

import sys
import itertools
import numpy as np

from triqs.lattice.tight_binding import TBLattice

from triqs.gf import Gf, MeshImFreq, Idx, inverse
from triqs.gf.gf_factories import make_gf_from_fourier

from triqs.operators import n, c, c_dag, Operator, dagger
from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.gw import get_gw_tensor

from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi_wr_from_chi_wk

from triqs_tprf.gw_solver import GWSolver

from gw_hubbard_dimer import GWHubbardDimer


class GWHubbardDimerMatrix:

    def __init__(self, beta=20.0, U=1.5, t=1.0, mu=0.0, nw=1024, maxiter=100):

        wmesh = MeshImFreq(beta, 'Fermion', nw)
            
        tb_opts = dict(
            units = [(1, 0, 0)],
            orbital_positions = [(0,0,0)] * 4,
            orbital_names = ['up_0', 'do_0', 'up_1', 'do_1'],
            )

        H0 = np.array([
            [ 0,  0,  1,  0],
            [ 0,  0,  0,  1],
            [ 1,  0,  0,  0],
            [ 0,  1,  0,  0],
            ])

        H_r = TBLattice(hopping = {
            (0,): -t * H0,
            }, **tb_opts)

        kmesh = H_r.get_kmesh(n_k=(1, 1, 1))
        self.e_k = H_r.fourier(kmesh)

        self.H_int = U*n('up',0) * n('do',0) + U*n('up',1) * n('do',1)

        self.fundamental_operators = \
            [c('up', 0), c('do', 0), c('up', 1), c('do', 1)]

        self.V_aaaa = get_gw_tensor(self.H_int, self.fundamental_operators)
        
        self.V_k = Gf(mesh=kmesh, target_shape=[4]*4)
        self.V_k.data[:] = self.V_aaaa    

        gw = GWSolver(self.e_k, self.V_k, wmesh, mu=mu)
        gw.solve_iter(maxiter=maxiter, gw=True, hartree=False, fock=False)
        gw.calc_real_space()
        
        self.gw = gw

        for key, val in gw.__dict__.items():
            setattr(self, key, val)


def test_gw_hubbard_dimer_matrix(verbose=False):

    beta = 20.0
    U = 1.5
    t = 1.0
    mu = 0.0
    
    gw_mat = GWHubbardDimerMatrix(
        beta = beta,
        U = U,
        t = t,
        mu = mu,
        nw = 1024,
        maxiter = 1,
        )

    gw = GWHubbardDimer(
        beta = beta,
        U = U,
        t = t,
        mu = mu,
        nw = 1024,
        maxiter = 1,
        )

    np.testing.assert_array_almost_equal(
        gw.P_wr[:, Idx(0, 0, 0)][0, 0, 0, 0].data,
        gw_mat.P_wr[:, Idx(0, 0, 0)][0, 0, 0, 0].data)

    np.testing.assert_array_almost_equal(
        gw.P_wr[:, Idx(1, 0, 0)][0, 0, 0, 0].data,
        gw_mat.P_wr[:, Idx(0, 0, 0)][0, 0, 2, 2].data)

    np.testing.assert_array_almost_equal(
        gw.W_wr[:, Idx(0, 0, 0)][0, 0, 0, 0].data,
        gw_mat.W_wr[:, Idx(0, 0, 0)][0, 0, 0, 0].data)

    np.testing.assert_array_almost_equal(
        gw.W_wr[:, Idx(1, 0, 0)][0, 0, 0, 0].data,
        gw_mat.W_wr[:, Idx(0, 0, 0)][0, 0, 2, 2].data)

    np.testing.assert_array_almost_equal(
        gw.sigma_wr[:, Idx(0, 0, 0)][0, 0].data,
        gw_mat.sigma_wr[:, Idx(0, 0, 0)][0, 0].data)

    np.testing.assert_array_almost_equal(
        gw.sigma_wr[:, Idx(1, 0, 0)][0, 0].data,
        gw_mat.sigma_wr[:, Idx(0, 0, 0)][0, 2].data)
    

if __name__ == '__main__':

    test_gw_hubbard_dimer_matrix()
