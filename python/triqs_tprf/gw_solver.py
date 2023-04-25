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

import numpy as np

import triqs.utility.mpi as mpi

from triqs.gf import Gf

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import gw_sigma, hartree_sigma, fock_sigma
from triqs_tprf.lattice import dynamical_screened_interaction_W
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk


class GWSolver():

    
    def __init__(self, e_k, V_k, wmesh, mu=None, g_wk=None):

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
        
        if g_wk is None:
            self.g_wk = self.g0_wk.copy()
        else:
            self.g_wk = g_wk
        
        self.sigma_wk = self.g0_wk.copy()
        self.sigma_wk.data[:] = 0.


    def calc_real_space(self):

        from triqs_tprf.lattice import fourier_wk_to_wr
        from triqs_tprf.lattice import chi_wr_from_chi_wk
        
        self.g0_wr = fourier_wk_to_wr(self.g0_wk)
        self.g_wr = fourier_wk_to_wr(self.g_wk)
        self.sigma_wr = fourier_wk_to_wr(self.sigma_wk)
        self.P_wr = chi_wr_from_chi_wk(self.P_wk)
        self.W_wr = chi_wr_from_chi_wk(self.W_wk)

        V_wk = self.W_wk.copy()
        bmesh = self.W_wk.mesh[0]
        V_wk.data[:] = 0
        for w in bmesh:
            V_wk[w, :] = self.V_k
        self.V_wk = V_wk
        self.V_wr = chi_wr_from_chi_wk(self.V_wk)
        
        
    def get_rho_loc(self, g_wk):
        
        wmesh = g_wk.mesh[0]
        kmesh = g_wk.mesh[1]
        g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
        g_w.data[:] = np.sum(g_wk.data, axis=1) / len(kmesh)
        rho = g_w.density()
        return rho
    
        
    def solve_iter(self, tol=1e-7, maxiter=100,
                   hartree=True, fock=True, gw=True,
                   verbose=True, spinless=False):

        if mpi.is_master_node():
            print(f'--> GWSolver.solve_iter')
            print(f'Hartree {hartree}')
            print(f'Fock    {fock}')
            print(f'GW      {gw}')
            print(f'spinless {spinless}')
            print()
        
        e_k, V_k, mu = self.e_k, self.V_k, self.mu

        g_wk = self.g_wk
        
        sigma_wk = g_wk.copy()
        sigma_wk.data[:] = 0.

        N_old = float('inf')
        g_wk_old = g_wk.copy()
        g_wk_old.data[:] = float('inf')
        
        for iter in range(maxiter):
            
            sigma_wk.data[:] = 0.

            if hartree:
                sigma_wk.data[:] += hartree_sigma(V_k, g_wk).data[None, ...]
            if fock:
                sigma_wk.data[:] += fock_sigma(V_k, g_wk).data[None, ...]

            if gw:
                P_wk = -imtime_bubble_chi0_wk(
                    g_wk, nw=len(self.wmesh)//2, verbose=False)

                if spinless:
                    W_wk = dynamical_screened_interaction_W(2*P_wk, V_k)
                else:
                    W_wk = dynamical_screened_interaction_W(P_wk, V_k)
                    
                W_dyn_wk, W_stat_k = split_into_dynamic_wk_and_constant_k(W_wk)
                
                np.testing.assert_array_almost_equal(W_stat_k.data, V_k.data)

                sigma_wk.data[:] += gw_sigma(W_dyn_wk, g_wk).data

            g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
            rho = self.get_rho_loc(g_wk)
            N = np.sum(np.diag(rho).real)

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
            self.W_dyn_wk = W_dyn_wk
            self.W_stat_k = W_stat_k
            self.P_wk = P_wk


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


def pade_analytical_continuation(
    g_wk, fmesh, n_points=32, freq_offset=0.05):
    
    g_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=g_wk.target_shape)
    wmesh = g_wk.mesh[0]

    for k in kmesh:
        
        if len(g_wk[:, k].data.shape) == 5:
            g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
            g_w.data[:] = g_wk[:, k].data[:, 0, 0, :, :]
        else:
            g_w = g_wk[:, k]

        g_f = g_fk[:, k]
        g_f.set_from_pade(g_w, n_points=n_points, freq_offset=freq_offset)
        
    return g_fk
