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

import itertools
import numpy as np
from scipy.linalg import block_diag

from triqs.gf import Gf, MeshImFreq
from triqs.operators import n, c, c_dag
from triqs.lattice.tight_binding import TBLattice

from triqs_tprf.OperatorUtils import operator_single_particle_transform
from triqs_tprf.OperatorUtils import quadratic_matrix_from_operator

from triqs_tprf.hf_solver import HartreeFockSolver

from triqs_tprf.gw import get_gw_tensor
from triqs_tprf.gw_solver import GWSolver


def test_gw_hubbard_dimer_hf():

    mu = 0.0
    beta = 1.0
    nw = 128

    t_u, t_d = 1.0, 0.5
    U, U_01 = 3.0, 1.0
    e_u0, e_d0 = 0.0, 1.0
    e_u1, e_d1 = 0.3, -0.3
    
    # Dimer with dens-dens interaction and hopping
    # vs. molecular orbitals with transformed interaction
    
    H_int = U * n('up',0) * n('do',0) + U * n('up',1) * n('do',1) + \
        U_01 * (n('up',0) + n('do',0)) * (n('up',1) + n('do',1))
    H_hop = \
        e_u0*n('up',0) + e_d0*n('do',0) + \
        e_u1*n('up',1) + e_d1*n('do',1) + \
        t_u * ( c_dag('up',0) * c('up',1) + c_dag('up',1) * c('up', 0) ) + \
        t_d * ( c_dag('do',0) * c('do',1) + c_dag('do',1) * c('do', 0) )
    H = H_hop + H_int
    
    print('-'*72)
    print(f'H_hop = {H_hop}')
    print(f'H_int = {H_int}')
    #print(f'H = {H}')

    # -- Single particle trasform to molecular orbitals
    
    U = 1 / np.sqrt(2) * np.array([
        [1,  1],
        [1, -1]])

    U = block_diag(U, U)
    np.testing.assert_array_almost_equal(U @ U.T.conj(), np.eye(4)) # Check unitarity
    
    gf_struct = [['up', 2], ['do', 2]]
    fundamental_operators = [c('up', 0), c('up',1), c('do', 0), c('do', 1)]
    
    Hp_hop = operator_single_particle_transform(H_hop, U, fundamental_operators)    
    Hp_int = operator_single_particle_transform(H_int, U, fundamental_operators)    
    Hp = Hp_hop + Hp_int

    print('-'*72)
    print(f'Hp_hop = {Hp_hop}')
    print(f'Hp_int = {Hp_int}')

    V_aaaa  = get_gw_tensor(H_int,  fundamental_operators)
    Vp_aaaa = get_gw_tensor(Hp_int, fundamental_operators)

    print('-'*72)
    print('--> V_aaaa')
    for i,j,k,l in itertools.product(range(4), repeat=4):
        value = V_aaaa[i, j, k, l]
        if np.abs(value) > 1e-9:
            print(f'{i}, {j}, {k}, {l} -- {value}')

    print('-'*72)
    print('--> Vp_aaaa')
    for i,j,k,l in itertools.product(range(4), repeat=4):
        value = Vp_aaaa[i, j, k, l]
        if np.abs(value) > 1e-9:
            print(f'{i}, {j}, {k}, {l} -- {value}')
            
    h_hop  = quadratic_matrix_from_operator(H_hop,  fundamental_operators)
    hp_hop = quadratic_matrix_from_operator(Hp_hop, fundamental_operators)

    print('-'*72)
    print(f'h_hop  =\n{h_hop.real}')
    print(f'hp_hop =\n{hp_hop.real}')
    np.testing.assert_array_almost_equal(U @ h_hop @ U.T.conj(), hp_hop) 
    
    # -- Lattices
    
    tb_opts = dict(
        units = [(1, 0, 0)],
        orbital_positions = [(0,0,0)] * 4,
        orbital_names = ['up_0', 'up_1', 'do_0', 'do_1'],
        )

    H_r = TBLattice(hopping = {(0,): h_hop,}, **tb_opts)
    Hp_r = TBLattice(hopping = {(0,): hp_hop,}, **tb_opts)

    kmesh = H_r.get_kmesh(n_k=(1, 1, 1))
    
    e_k = H_r.fourier(kmesh)
    ep_k = Hp_r.fourier(kmesh)

    # -- Hartree Fock solution
    
    hs  = HartreeFockSolver(e_k,  beta, H_int,  gf_struct)
    hsp = HartreeFockSolver(ep_k, beta, Hp_int, gf_struct)

    hs.solve_newton_mu(mu=mu)
    hsp.solve_newton_mu(mu=mu)
    print('-'*72)

    rho_aa  = hs.density_matrix()
    rhop_aa = hsp.density_matrix()

    np.testing.assert_array_almost_equal(
        hs.total_density(), hsp.total_density())
    
    np.testing.assert_array_almost_equal(U @ rho_aa @ U.T.conj(), rhop_aa) 

    # -- GW solution

    wmesh = MeshImFreq(beta, 'Fermion', nw)

    V_k = Gf(mesh=kmesh, target_shape=[4]*4)
    V_k.data[:] = V_aaaa    

    Vp_k = Gf(mesh=kmesh, target_shape=[4]*4)
    Vp_k.data[:] = Vp_aaaa    

    gw  = GWSolver(e_k,  V_k,  wmesh, mu=mu)
    gw.solve_iter(gw=False)

    gwp = GWSolver(ep_k, Vp_k, wmesh, mu=mu)
    gwp.solve_iter(gw=False)

    np.testing.assert_array_almost_equal(rho_aa, gw.rho)
    np.testing.assert_array_almost_equal(rhop_aa, gwp.rho)

    np.testing.assert_array_almost_equal(U @ gw.rho @ U.T.conj(), gwp.rho) 

    np.testing.assert_array_almost_equal(hs.M, np.squeeze(gw.sigma_wk.data[0, 0]))
    np.testing.assert_array_almost_equal(hsp.M, np.squeeze(gwp.sigma_wk.data[0, 0]))

    np.testing.assert_array_almost_equal(
        np.squeeze(gwp.sigma_wk.data[0, 0]),
        U @ np.squeeze(gw.sigma_wk.data[0, 0]) @ U.T.conj())



if __name__ == '__main__':

    test_gw_hubbard_dimer_hf()    
