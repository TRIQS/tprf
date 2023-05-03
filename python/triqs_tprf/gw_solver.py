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

from triqs.gf import Gf, MeshProduct, Idx
from triqs.gf.gf_factories import make_gf_from_fourier

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import rho_k_from_g_wk
from triqs_tprf.lattice import gw_dynamic_sigma, hartree_sigma, fock_sigma
from triqs_tprf.lattice import dynamical_screened_interaction_W
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr

from triqs_tprf.lattice import chi_wr_from_chi_wk
from triqs_tprf.lattice import chi_tr_from_chi_wr

from triqs_tprf.lattice import fourier_tr_to_wr
from triqs_tprf.lattice import fourier_wr_to_wk

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

from triqs_tprf.ParameterCollection import ParameterCollection

from triqs_tprf.ase_timing import Timer, timer


class GWSolver():

    
    def __init__(self, e_k, V_k, wmesh,
                 mu=None, g_wk=None, N_fix=None, N_tol=1e-5,
                 mu_bracket=None):

        self.timer = Timer()

        if mpi.is_master_node():
            print(self.logo())
            print(f'beta {wmesh.beta}')
            print(f'N_w {len(wmesh):6d}')
            print(f'N_k {len(e_k.mesh):6d}')
        
        self.e_k = e_k
        self.V_k = V_k
        self.V_r = make_gf_from_fourier(V_k)
        self.wmesh = wmesh
        self.N_fix = N_fix
        self.N_tol = N_tol

        self.mu = mu if mu is not None else 0.0
        
        if mu_bracket is None:
            self.mu_bracket = np.array([e_k.data.real.min(), e_k.data.real.max()])
        else:
            self.mu_bracket = mu_bracket

        self.g0_wk, self.mu = self.dyson_equation(mu, e_k, wmesh=wmesh, N_fix=N_fix)
        
        if g_wk is None:
            self.g_wk = self.g0_wk.copy()
        else:
            self.g_wk = g_wk
        
        self.sigma_wk = self.g0_wk.copy()
        self.sigma_wk.data[:] = 0.


    @timer('calc_real_space')
    def calc_real_space(self):

        if mpi.is_master_node():
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
        

    @timer('calc_real_freq')
    def calc_real_freq(self, fmesh, fbmesh=None, opts=dict()):

        if mpi.is_master_node():
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
            sigma_hartree_k = self.sigma_hartree_k,
            sigma_fock_k = self.sigma_fock_k,
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


    def calc_rho_loc(self, rho_r):
        rho_loc = np.array(rho_r[Idx(0, 0, 0)].data)
        return rho_loc

    
    def calc_total_density(self, rho_loc):
        N = np.sum(np.diag(rho_loc).real)
        return N
    
    
    @timer('Density matrix rho_r')
    def calc_rho_r(self, g_wk):
        rho_k = rho_k_from_g_wk(g_wk)
        rho_r = make_gf_from_fourier(rho_k)
        return rho_r

    
    @timer('Hartree Sigma_H')
    def hartree_sigma(self, V_k, rho_r):
        return make_gf_from_fourier(hartree_sigma(V_k, rho_r))

    
    @timer('Fock Sigma_F')
    def fock_sigma(self, V_r, rho_r):
        return make_gf_from_fourier(fock_sigma(V_r, rho_r))
    
        
    @timer('GW Sigma_dyn')
    def gw_dynamic_sigma(self, W_dyn_wk, g_wk):
 
        g_wr = fourier_wk_to_wr(g_wk)
        g_tr = fourier_wr_to_tr(g_wr)
        del g_wr

        W_dyn_wr = chi_wr_from_chi_wk(W_dyn_wk)
        W_dyn_tr = chi_tr_from_chi_wr(W_dyn_wr)
        del W_dyn_wr

        sigma_dyn_tr = gw_dynamic_sigma(W_dyn_tr, g_tr)
        del g_tr
        del W_dyn_tr

        sigma_dyn_wr = fourier_tr_to_wr(sigma_dyn_tr)
        del sigma_dyn_tr
        sigma_dyn_wk = fourier_wr_to_wk(sigma_dyn_wr)
        del sigma_dyn_wr

        return sigma_dyn_wk


    @timer('Polarization P_wk')
    def polarization(self, g_wk):
        P_wk = -imtime_bubble_chi0_wk(
            g_wk, nw=len(self.wmesh)//2, verbose=False)
        return P_wk


    @timer('Interaction W_wk')
    def screened_interaction(self, P_wk, V_k, spinless=False):
        if spinless:
            W_wk = dynamical_screened_interaction_W(2*P_wk, V_k)
        else:
            W_wk = dynamical_screened_interaction_W(P_wk, V_k)
        return W_wk

    
    @timer('Split W_dyn_k, W_stat_wk')
    def dynamic_and_static_interaction(self, W_wk):
        W_dyn_wk, W_stat_k = split_into_dynamic_wk_and_constant_k(W_wk)
        return W_dyn_wk, W_stat_k
    

    @timer('Dyson equation')
    def dyson_equation(self, mu, e_k, sigma_wk=None, wmesh=None, N_fix=None):

        if N_fix is None:
            g_wk = self._dyson_equation_dispatch(mu, e_k, sigma_wk=sigma_wk, wmesh=wmesh)
        else:
            # -- Seek chemical potential

            def target_function(mu):
                g_wk = self._dyson_equation_dispatch(mu, e_k, sigma_wk=sigma_wk, wmesh=wmesh)
                rho_r = self.calc_rho_r(g_wk)
                rho_loc = self.calc_rho_loc(rho_r)
                N = self.calc_total_density(rho_loc)
                return N - N_fix

            from scipy.optimize import root_scalar

            sol = root_scalar(target_function, method='brentq', bracket=self.mu_bracket, rtol=self.N_tol)
            mu = sol.root
            g_wk = self._dyson_equation_dispatch(mu, e_k, sigma_wk=sigma_wk, wmesh=wmesh)
            
        return g_wk, mu


    def _dyson_equation_dispatch(self, mu, e_k, sigma_wk=None, wmesh=None):

        if sigma_wk is not None:
            g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
        elif sigma_wk is None and wmesh is not None:
            g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
        else:
            raise NotImplementedError

        return g_wk
        
    
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
        else:
            verbose = False
        
        e_k, V_k, V_r, mu, N_fix = self.e_k, self.V_k, self.V_r, self.mu, self.N_fix

        g_wk = self.g_wk
        
        sigma_wk = g_wk.copy()
        sigma_wk.data[:] = 0.

        N_old = float('inf')
        g_wk_old = g_wk.copy()
        g_wk_old.data[:] = float('inf')
        
        for iter in range(maxiter):
            
            sigma_wk.data[:] = 0.

            if verbose: print('--> rho_r')
            rho_r = self.calc_rho_r(g_wk)

            if hartree:
                if verbose: print('--> Sigma Hartree')
                self.sigma_hartree_k = self.hartree_sigma(V_k, rho_r)
                sigma_wk.data[:] += self.sigma_hartree_k.data[None, ...]
            if fock:
                if verbose: print('--> Sigma Fock')
                self.sigma_fock_k = self.fock_sigma(V_r, rho_r)
                sigma_wk.data[:] += self.sigma_fock_k.data[None, ...]

            if gw:
                if verbose: print('--> Polarization')
                P_wk = self.polarization(g_wk)

                if verbose: print('--> Screened interaction')
                W_wk = self.screened_interaction(P_wk, V_k, spinless=spinless)
                    
                if verbose: print('--> Screened interaction split in static and dynamic')
                W_dyn_wk = W_wk.copy()
                for w in W_dyn_wk.mesh.components[0]:
                    W_dyn_wk[w,:] -= V_k

                if verbose: print('--> Sigma GW dynamic')
                self.sigma_dyn_wk = self.gw_dynamic_sigma(W_dyn_wk, g_wk)
                sigma_wk.data[:] += self.sigma_dyn_wk.data

            if verbose: print('--> Dyson equation')
            g_wk, mu = self.dyson_equation(mu, e_k, sigma_wk=sigma_wk, N_fix=N_fix)
            
            if verbose: print('--> done.')

            diff = np.max(np.abs(g_wk.data - g_wk_old.data))
            g_wk_old = g_wk

            if iter > 0 and verbose and mpi.is_master_node():
                print(f'GW: iter {iter:5d} max(abs(dg)) {diff:2.2E}')
                
            if diff < tol:
                break

        rho_r = self.calc_rho_r(g_wk)
        rho_loc = self.calc_rho_loc(rho_r)
        N = self.calc_total_density(rho_loc)

        self.rho_r = rho_r
        self.rho_loc = rho_loc
        self.N = N

        self.mu = mu
        self.g_wk = g_wk
        self.sigma_wk = sigma_wk

        if gw:
            self.W_wk = W_wk
            self.W_dyn_wk = W_dyn_wk
            self.P_wk = P_wk
            
        if mpi.is_master_node():
            print()
            self.timer.write()


    def get_local_density_matrix(self): return self.rho_loc

    def get_total_density(self): return self.N


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
