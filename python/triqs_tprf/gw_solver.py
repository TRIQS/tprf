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

from h5.formats import register_class 
import triqs.utility.mpi as mpi

from triqs.gf import Gf, MeshProduct
from triqs.gf.gf_factories import make_gf_from_fourier

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import rho_k_from_g_wk
from triqs_tprf.lattice import gw_sigma, hartree_sigma, fock_sigma
from triqs_tprf.lattice import dynamical_screened_interaction_W
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

from triqs_tprf.ParameterCollection import ParameterCollection

class GWSolver():

    
    def __init__(self, e_k, V_k, wmesh, mu=None, g_wk=None):

        if mpi.is_master_node():
            print(self.logo())
            print(f'beta {wmesh.beta}')
            print(f'N_w {len(wmesh):6d}')
            print(f'N_k {len(e_k.mesh):6d}')
        
        self.e_k = e_k
        self.V_k = V_k
        self.V_r = make_gf_from_fourier(V_k)
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

        print(f'--> GWSolver.calc_real_space')

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
        

    def calc_real_freq(self, fmesh, fbmesh=None, opts=dict()):

        print(f'--> GWSolver.calc_real_freq')

        self.g0_fk = pade_analytical_continuation(self.g0_wk, fmesh, **opts)
        self.g_fk = pade_analytical_continuation(self.g_wk, fmesh, **opts)
        self.sigma_fk = pade_analytical_continuation(self.sigma_wk, fmesh, **opts)

        gw_rf = ParameterCollection(
            mu = self.mu,
            e_k = self.e_k,
            V_k = self.V_k,
            g0_fk = self.g0_fk,
            g_fk = self.g_fk,
            )

        if fbmesh is None:
            fbmesh = fmesh
        
        if hasattr(self, 'P_wk'):
            self.P_fk = pade_analytical_continuation(self.P_wk, fbmesh, **opts)
            gw_rf.P_fk = self.P_fk
        if hasattr(self, 'W_wk'):
            self.W_fk = pade_analytical_continuation(self.W_wk, fbmesh, **opts)
            gw_rf.W_fk = self.W_fk
        if hasattr(self, 'W_dyn_wk'):
            self.W_dyn_fk = pade_analytical_continuation(self.W_dyn_wk, fbmesh, **opts)
            gw_rf.W_dyn_fk = self.W_dyn_fk

        return gw_rf


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
            print()
            print(f'--> GWSolver.solve_iter')
            print(f'spinless {spinless}')
            print(f' Hartree {hartree}')
            print(f'    Fock {fock}')
            print(f'      GW {gw}')
            print()
        
        e_k, V_k, V_r, mu = self.e_k, self.V_k, self.V_r, self.mu

        g_wk = self.g_wk
        
        sigma_wk = g_wk.copy()
        sigma_wk.data[:] = 0.

        N_old = float('inf')
        g_wk_old = g_wk.copy()
        g_wk_old.data[:] = float('inf')
        
        for iter in range(maxiter):
            
            sigma_wk.data[:] = 0.

            if verbose: print('--> rho_k')
            rho_k = rho_k_from_g_wk(g_wk)
            rho_r = make_gf_from_fourier(rho_k)

            if hartree:
                if verbose: print('--> Sigma Hartree')
                #self.sigma_hartree_k = hartree_sigma(V_k, g_wk)
                self.sigma_hartree_k = \
                    make_gf_from_fourier(hartree_sigma(V_r, rho_r))
                sigma_wk.data[:] += self.sigma_hartree_k.data[None, ...]
            if fock:
                if verbose: print('--> Sigma Fock')
                #self.sigma_fock_k = fock_sigma(V_k, g_wk)
                self.sigma_fock_k = \
                    make_gf_from_fourier(fock_sigma(V_r, rho_r))
                sigma_wk.data[:] += self.sigma_fock_k.data[None, ...]

            if gw:
                if verbose: print('--> Polarization')
                P_wk = -imtime_bubble_chi0_wk(
                    g_wk, nw=len(self.wmesh)//2, verbose=False)

                if verbose: print('--> Screened interaction')
                if spinless:
                    W_wk = dynamical_screened_interaction_W(2*P_wk, V_k)
                else:
                    W_wk = dynamical_screened_interaction_W(P_wk, V_k)
                    
                if verbose: print('--> Screened interaction split in static and dynamic')
                W_dyn_wk, W_stat_k = split_into_dynamic_wk_and_constant_k(W_wk)
                
                np.testing.assert_array_almost_equal(W_stat_k.data, V_k.data)

                if verbose: print('--> Sigma GW dynamic')
                self.sigma_dyn_wk = gw_sigma(W_dyn_wk, g_wk)
                sigma_wk.data[:] += self.sigma_dyn_wk.data

            if verbose: print('--> Dyson equation')
            g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
            rho = self.get_rho_loc(g_wk)
            N = np.sum(np.diag(rho).real)
            if verbose: print('--> done.')

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


    def __reduce_to_dict__(self):
        return self.__dict__

    
    @classmethod
    def __factory_from_dict__(cls, name, d):
        ret = cls(d['e_k'], d['V_k'], d['wmesh'])
        ret.__dict__.update(d)
        return ret

    
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


# -- Register ParameterCollection in Triqs formats
register_class(GWSolver)

    
def gf_tensor_to_matrix(g):
    
    N = g.target_shape[0] * g.target_shape[1] 
    M = g.target_shape[2] * g.target_shape[3] 
    g_mat = Gf(mesh=g.mesh, target_shape=[N, M])
    
    g_mat.data[:] = np.reshape(
        g.data, (g.data.shape[0], N, M))
    return g_mat


def gf_matrix_to_tensor(g_mat, target_shape):
    
    g = Gf(mesh=g_mat.mesh, target_shape=target_shape)
    shape = [g_mat.data.shape[0]] + list(target_shape)
    g.data[:] = np.reshape(g_mat.data, shape)
    return g


def pade_analytical_continuation(
    g_wk, fmesh, n_points=32, freq_offset=0.05):
    
    wmesh = g_wk.mesh[0]
    kmesh = g_wk.mesh[1]

    g_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=g_wk.target_shape)

    for k in kmesh:
        g_f = g_fk[:, k]
        g_w = g_wk[:, k]

        if len(g_wk.target_shape) == 4:
            g_w = gf_tensor_to_matrix(g_w)
            g_f = gf_tensor_to_matrix(g_f)
            
        g_f.set_from_pade(g_w, n_points=n_points, freq_offset=freq_offset)

        if len(g_wk.target_shape) == 4:
            g_f = gf_matrix_to_tensor(g_f, g_wk.target_shape)
            g_fk[:, k] = g_f
        
    return g_fk
