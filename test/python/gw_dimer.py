################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2032 by Hugo U.R. Strand
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
from scipy.linalg import block_diag

import triqs.utility.mpi as mpi

from triqs.lattice.tight_binding import TBLattice
from triqs.gf import Gf, MeshImFreq, Idx, MeshImTime
from triqs.operators import n, c, c_dag, Operator, dagger

from triqs_tprf.OperatorUtils import quadratic_matrix_from_operator
from triqs_tprf.OperatorUtils import operator_single_particle_transform

from triqs_tprf.hf_solver import HartreeSolver
from triqs_tprf.hf_solver import HartreeFockSolver

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import fourier_wk_to_wr

#from triqs_tprf.gw import bubble_PI_wk
#from triqs_tprf.gw import dynamical_screened_interaction_W
from triqs_tprf.lattice import gw_sigma, hartree_sigma, fock_sigma
#from triqs_tprf.gw import g0w_sigma
from triqs_tprf.gw import get_gw_tensor

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.gw import dynamical_screened_interaction_W_from_generalized_susceptibility
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k


def print_tensor(U, tol=1e-9):
    assert( len(U.shape) == 4)
    n = U.shape[0]
    
    for i,j,k,l in itertools.product(range(n), repeat=4):
        value = U[i, j, k, l]
        if np.abs(value) > tol:
            print(f'{i}, {j}, {k}, {l} -- {value}')


class GWSolver():

    
    def __init__(self, e_k, V_k, wmesh, mu=None):

        if mpi.is_master_node():
            print(self.logo())
            print(f'beta {wmesh.beta}')
            print(f'N_w {len(wmesh):6d}')
            print(f'N_k {len(e_k.mesh):6d}')
        
        self.e_k = e_k
        self.V_k = V_k
        self.wmesh = wmesh

        self.mu = mu if mu is not None else 0.0

        self.g0_wk = lattice_dyson_g0_wk(mu=self.mu, e_k=e_k, mesh=wmesh)
        
        self.g_wk = self.g0_wk.copy()
        self.sigma_wk = self.g0_wk.copy()
        self.sigma_wk.data[:] = 0.


    def get_rho(self, g_wk):
        
        wmesh = g_wk.mesh[0]
        kmesh = g_wk.mesh[1]
        g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
        g_w.data[:] = np.sum(g_wk.data, axis=1) / len(kmesh)
        rho = g_w.density()
        return rho
    
        
    def solve_iter(self, tol=1e-7, maxiter=100,
                   hartree=True, fock=True, gw=True,
                   verbose=True):

        e_k, V_k, mu = self.e_k, self.V_k, self.mu

        g_wk = self.g_wk
        
        sigma_wk = g_wk.copy()
        sigma_wk.data[:] = 0.

        N_old = float('inf')
        g_wk_old = g_wk.copy()
        g_wk_old.data[:] = float('inf')
        
        for iter in range(maxiter):

            if gw:
                chi00_wk = imtime_bubble_chi0_wk(g_wk, nw=len(self.wmesh)//2)

                W_wk = \
                    dynamical_screened_interaction_W_from_generalized_susceptibility(
                        -chi00_wk, V_k)
                W_dyn_wk, W_stat_k = split_into_dynamic_wk_and_constant_k(W_wk)
            
            sigma_wk.data[:] = 0.
            
            if hartree:
                sigma_wk.data[:] += hartree_sigma(V_k, g_wk).data[None, ...]
            if fock:
                sigma_wk.data[:] += fock_sigma(V_k, g_wk).data[None, ...]

            if gw:
                sigma_wk.data[:] += gw_sigma(W_wk, g_wk).data

            g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
            rho = self.get_rho(g_wk)
            N = np.sum(np.diag(rho).real)

            #dN = np.abs(N_old - N)
            #N_old = N
            diff = np.max(np.abs(g_wk.data - g_wk_old.data))
            g_wk_old = g_wk

            if iter > 0 and verbose and mpi.is_master_node():
                print(f'GW: iter {iter:5d} max(abs(dg)) {diff:2.2E}')
                
            if diff < tol:
                break

        self.N = N
        self.rho = rho

        self.g_wk = g_wk
        self.sigma_wk = sigma_wk

        if gw:
            self.W_wk = W_wk
            self.chi00_wk = chi00_wk


    def logo(self):
        if 'UTF' in sys.stdout.encoding:
            logo = """
╔╦╗╦═╗╦╔═╗ ╔═╗  
 ║ ╠╦╝║║═╬╗╚═╗  
 ╩ ╩╚═╩╚═╝╚╚═╝  GW
TRIQS: GW solver
"""
        else:
            logo = """
 _____ ___ ___ ___  ___
|_   _| _ \_ _/ _ \/ __| 
  | | |   /| | (_) \__ \ 
  |_| |_|_\___\__\_\___/  GW

    TRIQS: GW solver
"""
        return logo


def hubbard_atom_hf():

    beta = 1.0
    mu = 0.0
    nw = 128
    
    U = 0.1
    t = 1.0
    e_u, e_d = 0.0, 0.0

    H_int = U * n('up',0) * n('do',0)
    H_hop = t * (c_dag('up',0) * c('do',0) + c_dag('do',0) * c('up',0)) + \
        e_u * n('up', 0) + e_d * n('do', 0)
    
    print(f'H_hop =\n{H_hop}')
    print(f'H_int =\n{H_int}')

    U = 1 / np.sqrt(2) * np.array([
        [1,  1],
        [1, -1]])

    np.testing.assert_array_almost_equal(U @ U.T.conj(), np.eye(2)) # Check unitarity
    
    gf_struct = [['up', 1], ['do', 1]]
    fundamental_operators = [c('up', 0), c('do', 0)]
    
    Hp_hop = operator_single_particle_transform(H_hop, U, fundamental_operators)    
    Hp_int = operator_single_particle_transform(H_int, U, fundamental_operators)    
    Hp = Hp_hop + Hp_int

    print('-'*72)
    print(f'Hp_hop = {Hp_hop}')
    print(f'Hp_int = {Hp_int}')
    
    h_hop  = quadratic_matrix_from_operator(H_hop,  fundamental_operators)
    hp_hop = quadratic_matrix_from_operator(Hp_hop, fundamental_operators)

    print('-'*72)
    print(f'h_hop  =\n{h_hop.real}')
    print(f'hp_hop =\n{hp_hop.real}')
    np.testing.assert_array_almost_equal(U @ h_hop @ U.T.conj(), hp_hop) 
    
    tb_opts = dict(
        units = [(1, 0, 0)],
        orbital_positions = [(0,0,0)] * 2,
        orbital_names = ['up_0', 'do_0'],
        )

    H_r = TBLattice(hopping = {(0,): h_hop,}, **tb_opts)
    Hp_r = TBLattice(hopping = {(0,): hp_hop,}, **tb_opts)

    kmesh = H_r.get_kmesh(n_k=(1, 1, 1))
    
    e_k = H_r.fourier(kmesh)
    ep_k = Hp_r.fourier(kmesh)

    # -- Hartree Fock? solution
    
    hs  = HartreeFockSolver(e_k,  beta, H_int,  gf_struct)
    hsp = HartreeFockSolver(ep_k, beta, Hp_int, gf_struct)

    hs.solve_newton_mu(mu=mu)
    hsp.solve_newton_mu(mu=mu)
    print('-'*72)

    rho_aa  = hs.density_matrix()
    rhop_aa = hsp.density_matrix()

    np.testing.assert_array_almost_equal(U @ rho_aa @ U.T.conj(), rhop_aa)
    
    np.testing.assert_array_almost_equal(
        hs.total_density(), hsp.total_density())

    # -- GW solution

    V_aaaa = get_gw_tensor(H_int, fundamental_operators)
    Vp_aaaa = get_gw_tensor(Hp_int, fundamental_operators)

    print('-'*72)
    print('--> V_aaaa')
    print_tensor(V_aaaa)
    print('--> Vp_aaaa')
    print_tensor(Vp_aaaa)

    V_k = Gf(mesh=kmesh, target_shape=[2]*4)
    V_k.data[:] = V_aaaa    

    Vp_k = Gf(mesh=kmesh, target_shape=[2]*4)
    Vp_k.data[:] = Vp_aaaa    

    wmesh = MeshImFreq(beta, 'Fermion', nw)
    
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
    
            
def hubbard_dimer_hf():

    mu = 0.0
    beta = 1.0
    nw = 128

    if True:
        t_u, t_d = 1.0, 0.5
        U, U_01 = 3.0, 1.0
        e_u0, e_d0 = 0.0, 1.0
        e_u1, e_d1 = 0.3, -0.3
    else:
        t_u, t_d = 1.0, 1.0
        U, U_01 = 1.0, 1.0
        e_u0, e_d0 = 0.0, 0.0
        e_u1, e_d1 = 0.0, 0.0
    
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

    #exit()
            
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

    # -- Hartree Fock? solution
    
    hs  = HartreeFockSolver(e_k,  beta, H_int,  gf_struct)
    hsp = HartreeFockSolver(ep_k, beta, Hp_int, gf_struct)

    hs.solve_newton_mu(mu=mu)
    hsp.solve_newton_mu(mu=mu)
    print('-'*72)

    rho_aa  = hs.density_matrix()
    rhop_aa = hsp.density_matrix()

    np.testing.assert_array_almost_equal(
        hs.total_density(), hsp.total_density())

    #exit()
    
    np.testing.assert_array_almost_equal(U @ rho_aa @ U.T.conj(), rhop_aa) 

    # -- GW solution

    wmesh = MeshImFreq(beta, 'Fermion', nw)
    g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)    
    gp_wk = lattice_dyson_g0_wk(mu=mu, e_k=ep_k, mesh=wmesh)    

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
    

def hubbard_dimer_gw():

    beta = 10.0
    U = 1.0
    t = 1.0
    nw = 1024
    mu = 0.0
    
    tb_opts = dict(
        units = [(1, 0, 0)],
        orbital_positions = [(0,0,0)] * 2,
        orbital_names = ['up_0', 'do_0'],
        )

    H_r = TBLattice(hopping = {
        (+1,): -t * np.eye(2),
        (-1,): -t * np.eye(2),
        }, **tb_opts)

    kmesh = H_r.get_kmesh(n_k=(2, 1, 1))
    e_k = H_r.fourier(kmesh)
    print(e_k.data)

    H_int = U * n('up',0) * n('do',0)
    print(f'H_int = {H_int}')
    
    fundamental_operators = [c('up', 0), c('do', 0)]
    V_aaaa = get_gw_tensor(H_int, fundamental_operators)

    print('V_aaaa =')
    print_tensor(V_aaaa)

    V_k = Gf(mesh=kmesh, target_shape=[2]*4)
    V_k.data[:] = V_aaaa    

    wmesh = MeshImFreq(beta, 'Fermion', nw)
    gw = GWSolver(e_k, V_k, wmesh, mu=mu)
    g0_wk = gw.g_wk.copy()
    g0_wr = fourier_wk_to_wr(g0_wk)

    if False:
        #from triqs_tprf.lattice import fourier_wk_to_wr
        from triqs_tprf.lattice import fourier_wr_to_tr
        from triqs_tprf.lattice import chi0_tr_from_grt_PH
        from triqs_tprf.lattice import chi_wr_from_chi_tr
        from triqs_tprf.lattice import chi_wk_from_chi_wr

        g_wk = g0_wk
        nw = len(g_wk.mesh.components[0]) // 2    
        g_wr = fourier_wk_to_wr(g_wk)
        g_tr = fourier_wr_to_tr(g_wr)
        chi00_tr = chi0_tr_from_grt_PH(g_tr)
        PI_tr = -chi00_tr # NB! sign difference
        PI_wr = chi_wr_from_chi_tr(PI_tr, nw=nw)
        PI_wk = chi_wk_from_chi_wr(PI_wr)

    
    gw.solve_iter(maxiter=1, gw=True, hartree=False, fock=False)
    g_wk = gw.g_wk

    from triqs_tprf.lattice import chi_wr_from_chi_wk
    P_wr = -chi_wr_from_chi_wk(gw.chi00_wk)
    W_wr = chi_wr_from_chi_wk(gw.W_wk)

    #np.testing.assert_array_almost_equal(PI_wk.data, -gw.chi00_wk.data)
    
    # Non-interacting reference Gf

    # Pina Romaniello
    # Hubbard Dimer in GW and Beyond
    # https://www.cond-mat.de/events/correl21/manuscripts/correl21.pdf    

    # E. Pavarini and E. Koch (eds.)
    # Simulating Correlations with Computers
    # Modeling and Simulation Vol. 11
    # Forschungszentrum Ju ̈lich, 2021, ISBN 978-3-95806-529-1
    # http://www.cond-mat.de/events/correl21
        
    g_0_w = Gf(mesh=wmesh, target_shape=[])
    g_1_w = Gf(mesh=wmesh, target_shape=[])

    # Eq. 18
    
    for w in wmesh:
        g_0_w[w] = +0.5/(w - 2*t) + 0.5/(w + 2*t)
        g_1_w[w] = -0.5/(w - 2*t) + 0.5/(w + 2*t)

    np.testing.assert_array_almost_equal(
        g0_wr[:, Idx(0, 0, 0)][0, 0].data, g_0_w.data)

    np.testing.assert_array_almost_equal(
        g0_wr[:, Idx(1, 0, 0)][0, 0].data, g_1_w.data)

    # Eq. 22

    bmesh = P_wr.mesh[0]

    P_0_w = Gf(mesh=bmesh, target_shape=[])
    P_1_w = Gf(mesh=bmesh, target_shape=[])
    
    for w in bmesh:
        P_0_w[w] = + 0.25 / (w - 4*t) - 0.25 / (w + 4*t)
        P_1_w[w] = - 0.25 / (w - 4*t) + 0.25 / (w + 4*t)

    np.testing.assert_array_almost_equal(
        P_wr[:, Idx(0, 0, 0)][0, 0, 0, 0].data, P_0_w.data)

    np.testing.assert_array_almost_equal(
        P_wr[:, Idx(1, 0, 0)][0, 0, 0, 0].data, P_1_w.data)
        
    W_0_w = Gf(mesh=bmesh, target_shape=(1,))
    W_1_w = Gf(mesh=bmesh, target_shape=(1,))

    h2 = 4*(2*t)**2 + 2*U*(2*t)
    
    for w in bmesh:
        W_0_w[w] = U + U**2 * 2*t/(complex(w)**2 - h2)
        W_1_w[w] = - U**2 * 2*t /(complex(w)**2 - h2)
        
    from triqs.plot.mpl_interface import oplot, oploti, oplotr, plt

    plt.figure(figsize=(9, 9))

    subp = [3, 2, 1]
    xlim = [-10, 10]
    
    plt.subplot(*subp); subp[-1] += 1
    oplot(g0_wr[:, Idx(0, 0, 0)][0, 0], 'k-x')
    oplot(g_0_w, 'g-+')
    plt.xlim(xlim)
    
    plt.subplot(*subp); subp[-1] += 1
    oplot(g0_wr[:, Idx(1, 0, 0)][0, 0], 'k-x')
    oplot(g_1_w, 'g-+')
    plt.xlim(xlim)

    plt.subplot(*subp); subp[-1] += 1
    oplot(P_wr[:, Idx(0, 0, 0)][0, 0, 0, 0], 'k-x')
    oplot(P_0_w, 'g-+')
    plt.xlim(xlim)

    plt.subplot(*subp); subp[-1] += 1
    oplot(P_wr[:, Idx(1, 0, 0)][0, 0, 0, 0], 'k-x')
    oplot(P_1_w, 'g-+')
    plt.xlim(xlim)

    plt.subplot(*subp); subp[-1] += 1
    oplot(W_wr[:, Idx(0, 0, 0)][0, 0, 0, 0], 'k-x')
    oplot(W_0_w - U, 'g-+')
    plt.xlim(xlim)

    plt.subplot(*subp); subp[-1] += 1
    oplot(W_wr[:, Idx(1, 0, 0)][0, 0, 0, 0], 'k-x')
    oplot(W_1_w, 'g-+')
    plt.xlim(xlim)

    plt.tight_layout()
    plt.show()
    
if __name__ == '__main__':

    #hubbard_atom_hf()
    #hubbard_dimer_hf()
    hubbard_dimer_gw()
