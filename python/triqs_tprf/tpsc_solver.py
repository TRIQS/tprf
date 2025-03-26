# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2025 by Xaver Landerl, H. U.R. Strand
# Authors: Xaver Landerl, H. U.R. Strand
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

import numpy as np

from scipy.optimize import brentq

from triqs.gf import Gf, MeshProduct
from triqs.gf import MeshCycLat, MeshBrZone
from triqs.gf import MeshImFreq, MeshDLRImFreq
from triqs.gf import MeshImTime, MeshDLRImTime
from triqs.gf import density as density_from_gf

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import gw_dynamic_sigma

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import fourier_wr_to_tr
from triqs_tprf.lattice import fourier_tr_to_wr
from triqs_tprf.lattice import fourier_wr_to_wk

import sys

def tpsc_banner():
    if 'UTF' in sys.stdout.encoding:
        # https://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20tpsc
        banner = r"""
╔╦╗╦═╗╦╔═╗ ╔═╗  ┌┬┐┌─┐┌─┐┌─┐
 ║ ╠╦╝║║═╬╗╚═╗   │ ├─┘└─┐│  
 ╩ ╩╚═╩╚═╝╚╚═╝   ┴ ┴  └─┘└─┘
Two-Particle Self-Consistent method"""
    else:
        # https://patorjk.com/software/taag/#p=display&f=Small&t=TRIQS%20TPSC
        banner = r"""
  _____ ___ ___ ___  ___   _____ ___  ___  ___ 
 |_   _| _ \_ _/ _ \/ __| |_   _| _ \/ __|/ __|
   | | |   /| | (_) \__ \   | | |  _/\__ \ (__ 
   |_| |_|_\___\__\_\___/   |_| |_|  |___/\___|

  Two-Particle Self-Consistent method"""
    return banner



def fourier_wk_to_tr(g_wk):
    """
    Fourier-transforms a Green's function from wk to tr representation.

    Requires
    --------
    g_wk must be defined on a MeshProduct(Mesh(DLR)ImFreq, MeshBrZone)

    Parameters
    ----------
    triqs.gf g_wk       : TRIQS Green's function

    Returns
    -------
    triqs.gf g_tr       : TRIQS Green's function on a MeshProduct(Mesh(DLR)ImTime, MeshCycLat)
    """

    # extract mesh
    mesh = g_wk.mesh

    # check input
    if not isinstance(mesh, MeshProduct):
        raise TypeError('g_wk.mesh must be of type \'MeshProduct\'!')
    
    # extract meshes
    w_mesh, k_mesh = mesh.components

    if not isinstance(w_mesh, MeshDLRImFreq) or isinstance(w_mesh, MeshImFreq):
        raise TypeError('g_wk.mesh.components[0] must be of type \'Mesh(DLR)ImFreq\'!')
    
    if not isinstance(k_mesh, MeshBrZone):
        raise TypeError('g_wk.mesh.components[1] must be of type \'MeshBrZone\'!')
    
    # fourier transform
    g_wr = fourier_wk_to_wr(g_wk)
    g_tr = fourier_wr_to_tr(g_wr)

    # return result
    return g_tr


def fourier_wk_to_mtr(g_wk):
    """
    Fourier-transforms a Green's function from wk to (-t,-r) representation.

    Requires
    --------
    g_wk must be defined on a MeshProduct(Mesh(DLR)ImFreq, MeshBrZone)

    Parameters
    ----------
    triqs.gf g_wk       : TRIQS Green's function

    Returns
    -------
    triqs.gf g_mtr      : TRIQS Green's function on a MeshProduct(Mesh(DLR)ImTime, MeshCycLat)
    """

    # extract mesh
    mesh = g_wk.mesh

    # check input
    if not isinstance(mesh, MeshProduct):
        raise TypeError('g_wk.mesh must be of type \'MeshProduct\'!')
    
    # extract meshes
    w_mesh, k_mesh = mesh.components

    if not isinstance(w_mesh, MeshDLRImFreq) or isinstance(w_mesh, MeshImFreq):
        raise TypeError('g_wk.mesh.components[0] must be of type \'Mesh(DLR)ImFreq\'!')
    
    if not isinstance(k_mesh, MeshBrZone):
        raise TypeError('g_wk.mesh.components[1] must be of type \'MeshBrZone\'!')
    
    # fourier transform
    g_wk_conj = g_wk.conjugate()
    g_wr_conj = fourier_wk_to_wr(g_wk_conj)
    g_tr_conj = fourier_wr_to_tr(g_wr_conj)
    g_mtr = g_tr_conj.conjugate()

    # return result
    return g_mtr


def fourier_tr_to_wk(g_tr):
    """
    Fourier-transforms a Green's function from wk to tr representation.

    Requires
    --------
    g_tr must be defined on a MeshProduct(Mesh(DLR)ImTime, MeshCycLat)

    Parameters
    ----------
    triqs.gf g_tr       : TRIQS Green's function

    Returns
    -------
    triqs.gf g_wk       : TRIQS Green's function on a MeshProduct(Mesh(DLR)ImFreq, MeshBrZone)
    """

    # extract mesh
    mesh = g_tr.mesh

    # check input
    if not isinstance(mesh, MeshProduct):
        raise TypeError('g_tr.mesh must be of type \'MeshProduct\'!')
    
    # extract meshes
    w_mesh, k_mesh = mesh.components

    if not isinstance(w_mesh, MeshDLRImTime) or isinstance(w_mesh, MeshImTime):
        raise TypeError('g_tr.mesh.components[0] must be of type \'Mesh(DLR)ImTime\'!')
    
    if not isinstance(k_mesh, MeshCycLat):
        raise TypeError('g_tr.mesh.components[1] must be of type \'MeshCycLat\'!')
    
    # fourier transform
    g_wr = fourier_tr_to_wr(g_tr)
    g_wk = fourier_wr_to_wk(g_wr)

    # return result
    return g_wk



class tpsc_solver:
    """
    Class to calculate the TPSC-approximation to the single-band Hubbard Model, introduced in [1].

    Can calculate TPSC-susceptibilities, self-energy and Green's function.

    [1] Non-perturbative many-body approach to the Hubbard model and single-particle pseudogap, Y.M. Vilk & A.-M.S. Tremblay
    """

    def __init__(self, n, U, wmesh, e_k, docc='tpsc_ansatz', verbose=True):
        """
        Initialize a tpsc_solver object.

        Parameters
        ----------
        n : double
            total electron density
        U : double
            Hubbard interaction
        wmesh : MeshImFreq or MeshDLRImFreq
                imaginary frequency mesh
        e_k : Gf
              dispersion relation
        docc : double, optional
               double occupancy (default: use TPSC-Ansatz)
        verbose : bool, optional
                  if true, print out information (default: True)
        """
        
        self.n = n
        self.U = U
        self.wmesh = wmesh
        self.beta = self.wmesh.beta
        self.e_k = e_k

        if docc == 'tpsc_ansatz':
            self.use_tpsc_ansatz = True
        else:
            assert( type(docc) == float )
            self.use_tpsc_ansatz = False
            self.docc = docc

        self.verbose = verbose

        self.vprint(tpsc_banner()+'\n')
        if self.use_tpsc_ansatz == True:
            self.vprint("  Using the TPSC-Ansatz for the double occupancy.")
        else:
            self.vprint(f'docc = {self.docc} (not using the TPSC-ansatz)')
        self.vprint()
        self.vprint(f'len(wmesh) = {len(wmesh)}')
        self.vprint(f'nk = {len(self.e_k.mesh)}')
        self.vprint(f'beta = {self.beta}')
        self.vprint(f'n = {self.n}')
        self.vprint(f'U = {self.U}')
    
    def vprint(self, string=None):
        if self.verbose == True:
            if string is not None:
                print(string)
            else:
                print()

    def solve(self, calc_sigma=True, calc_g=True, check_self_consistency=True,
              Usp_tol=None, Uch_tol=None, Uch_max=100.):
        """
        Runs the TPSC-Calculation on the specified model.

        Parameters
        ----------
        calc_sigma : bool, optional
                     if True, the TPSC-self-energy is calculated
                     (default: True)
        calc_g : bool, optional
                 if True, the lattice Green's function is calculated
                 (default: True)
        check_self_consistency : bool, optional
                                 if True, self-consistency of the result is checked
                                 (default: True)
        Usp_tol : double, optional
                  the computed Usp will satisfy np.allclose(Usp, Usp_exact, atol=Usp_tol)
                  (default: 2e-12, same as scipy.optimize.brentq)
        Uch_tol : double, optional
                  the computed Uch will satisfy np.allclose(Uch, Uch_exact, atol=Usp_tol)
                  (default: 2e-12, same as scipy.optimize.brentq)
        Uch_max : maximum value in Uch search
                  (default: 100.)
        """

        # check consistency of input
        if calc_g == True:
            if calc_sigma != True:
                raise ValueError("Sigma must known to calculate the Green's function!")
        if check_self_consistency == True:
            if calc_sigma != True or calc_g != True:
                raise ValueError("Sigma and G must be known to check self-consistency!")
              
        self.vprint("\nCalculating non-interacting quantities...")
        
        self.g0_wk = self._calc_g0_wk(self.n)

        # clean up tprf api FIXME
        if isinstance(self.wmesh, MeshDLRImFreq):
            self.chi0_wk = 2*imtime_bubble_chi0_wk(self.g0_wk, nw=-1, verbose=False)

        elif isinstance(self.wmesh, MeshImFreq):
            self.chi0_wk = 2*imtime_bubble_chi0_wk(self.g0_wk, nw=self.wmesh.n_iw, verbose=False)

        self.vprint("\nCalculating first level of approximation...")
        self._calc_first_level_approx(Usp_tol=Usp_tol, Uch_tol=Uch_tol, Uch_max=Uch_max)
        
        if calc_sigma == True:
            self.vprint()
            self.vprint("Calculating self-energy...")
            self._get_sigma()
        
        if calc_g == True:
            self.vprint()
            self.vprint("Calculating Green\'s function...")
            self.g_wk, self.mu = self._calc_g_wk(self.n, self.sigma_wk)

        if check_self_consistency == True:
            self.vprint()
            self.vprint("Checking self-consistency requirements...")
            self._check_for_self_consistency()
        
        self.vprint()
        self.vprint(" DONE! ")

    def _calc_first_level_approx(self, Usp_tol=None, Uch_tol=None, Uch_max=100.):
        """
        Runs the first-level approximation of the TPSC-Calculation on the specified model.
        If docc is specified in the model, it will be used to evaluate the sum rules.
        If docc is not specified in the model, the TPSC-Ansatz will be used to evaluate the sum rules.

        Parameters
        ----------
        Usp_tol : double, optional
                  the computed Usp will satisfy np.allclose(Usp, Usp_exact, atol=Usp_tol)
                  (default: 2e-12, same as scipy.optimize.brentq)
        Uch_tol : double, optional
                  the computed Uch will satisfy np.allclose(Uch, Uch_exact, atol=Usp_tol)
                  (default: 2e-12, same as scipy.optimize.brentq)
                  (if Usp_tol is set and Uch_tol is not, then Uch_tol=Usp_tol)
        Uch_max : maximum value in Uch search
                  (default: 100.)
        """
        
        # calculate vertices

        # calculate Usp
        self.vprint()
        self.vprint("   Calculating Usp...")
        def _Usp_root(Usp):

            chi_wk = self._solve_rpa(Usp)
            LHS = self._get_density(chi_wk)

            if self.use_tpsc_ansatz == True:
                RHS = self.n - Usp/self.U*self.n*self.n/2
            else:
                RHS = self.n - 2*self.docc

            return LHS - RHS
        
        # here chi_sp diverges (1e-7 for numerical stability)
        Usp_max = 2.0/np.amax(self.chi0_wk.data).real - 1e-7
        self.Usp = brentq(_Usp_root, 0.0, Usp_max, **({'xtol':Usp_tol} if Usp_tol is not None else {}))
        
        # calculate Uch
        self.vprint()
        self.vprint("   Calculating Uch...")
        def _Uch_root(Uch):
            
            chi_wk = self._solve_rpa(-Uch)
            LHS = self._get_density(chi_wk)
            
            if self.use_tpsc_ansatz == True:
                RHS = self.n + self.Usp/self.U*self.n*self.n/2 - self.n*self.n
            else:
                RHS = self.n + 2*self.docc - self.n*self.n
            
            return LHS - RHS
        
        self.Uch = brentq(_Uch_root, 0.0, Uch_max, **({'xtol':Uch_tol} if Uch_tol is not None else {}))

        # calculate susceptibilities
        self.vprint()
        self.vprint("   Calculating chisp_wk and chich_wk...")
        self.chisp_wk = self._solve_rpa(self.Usp)
        self.chich_wk = self._solve_rpa(-self.Uch)
        
        # calculate double occupation
        if self.use_tpsc_ansatz == True:
            self.vprint()
            self.vprint("   Calculating double occupancy...")
            self.docc = self.Usp/self.U*self.n*self.n/4
        
        # print out results
        self.vprint()
        self.vprint("Summary first level approximation:")
        self.vprint("    Usp = " + str(self.Usp) + ", Uch = " + str(self.Uch))
        self.vprint("    Double Occupation <n_up*n_down> = " + str(self.docc))

    def _check_for_self_consistency(self):
        """
        Checks the evaluated model for self-consistency.
        
        Requires
        --------
        Must have calculated chi_sp and chi_ch.
        
        Parameters
        ----------
        self    : self

        Returns
        -------
        None

        """
        
        # check that the sum rule is fulfilled
        self.vprint()
        self.vprint("------------------------------------------------------------------------------------------")
        self.vprint('Doing self-consistency check of first-level approximation...')
        check_sum_rule = self._get_density(self.chisp_wk + self.chich_wk) - (2*self.n - self.n**2)
        self.vprint(f'The sum rule is fulfilled with an accuracy of {abs(check_sum_rule)}.')
        self.vprint("------------------------------------------------------------------------------------------")
        self.vprint('Doing self-consistency check of second-level approximation...')

        # get traces (the second term comes from the HF self-energy)
        trace_F1 = self._get_density(self.sigma_wk * self.g0_wk) / self.U
        trace_F2 = self._get_density(self.sigma_wk * self.g_wk) / self.U

        # get the relative difference
        rel_diff = np.abs((trace_F1 - trace_F2) / trace_F1) * 100

        # check consistency
        self.vprint('Sum_k {Sigma^(2)(k) * G^(1)(k)} - U*<n_up * n_down> = ' + str(np.abs(trace_F1 - self.docc)))
        self.vprint('Sum_k {Sigma^(2)(k) * G^(2)(k)} - U*<n_up * n_down> = ' + str(float(np.abs(trace_F2 - self.docc))))
        self.vprint(f'Relative difference of traces = {rel_diff:.2f}%')

        
    def _solve_rpa(self, U_vert):
        """
        Calculates an rpa-susceptibility with a scalar vertex from the non-interacting susceptibility.

        Parameters
        U_vert : double
                 RPA vertex

        Returns
        -------
        
        """

        ### FIXME
        if False:
            V = 0.5 * U_vert * np.ones((1, 1, 1, 1), dtype=complex)
            chi_wk = solve_rpa_PH(self.chi0_wk, V)
            return chi_wk

        chi_wk = self.chi0_wk.copy()

        # fill with data
        chi_wk.data[:] = self.chi0_wk.data[:]/(1 - U_vert/2*self.chi0_wk.data[:])

        # return results
        return chi_wk

    
    def _get_density(self, g_wk):

        wmesh = g_wk.mesh[0]
        nk = g_wk.data.shape[1]
        g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
        g_w.data[:] = np.sum(g_wk.data, axis=1) / nk

        assert( g_w.target_shape == tuple([1]*len(g_w.target_shape)) )
    
        g_w_scalar = Gf(mesh=g_w.mesh, target_shape=[])
        g_w_scalar.data[:] = np.squeeze(g_w.data)
        
        dens = density_from_gf(g_w_scalar).real
        
        if g_w.mesh.statistic == 'Boson':
            dens = -dens

        return dens   

    def _calc_g0_wk(self, target_density):

        # find the mu that leads to the correct density
        mu_min, mu_max = np.min(self.e_k.data.real), np.max(self.e_k.data.real)

        def target_function(mu):
            g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=self.e_k, mesh=self.wmesh)
            density = self._get_density(g0_wk)
            return target_density/2 - density
    
        mu0 = brentq(target_function, mu_min, mu_max)
        g0_wk = lattice_dyson_g0_wk(mu=mu0, e_k=self.e_k, mesh=self.wmesh)
        return g0_wk


    def _calc_g_wk(self, target_density, sigma_wk):

        # find the mu that leads to the correct density
        mu_min, mu_max = np.min(self.e_k.data.real), np.max(self.e_k.data.real)

        def target_function(mu):
            g_wk = lattice_dyson_g_wk(mu, self.e_k, sigma_wk)
            density = self._get_density(g_wk)
            return target_density/2 - density
    
        mu = brentq(target_function, mu_min, mu_max)
        g_wk = lattice_dyson_g_wk(mu, self.e_k, sigma_wk)
        
        return g_wk, mu
    
    
    def _calc_sigma_deprecated(self):
        """
        Calculates the second-level approximation of the self-energy.
        """

        # define effective potential
        V_wk = self.U/8*(3*self.Usp*self.chisp_wk + self.Uch*self.chich_wk)

        # get V(-t,-r)
        V_mtr = fourier_wk_to_mtr(V_wk)

        # get G(t,r)
        g0_tr = fourier_wk_to_tr(self.g0_wk)

        # multiply V(-t,-r) * G0(t,r) = Sigma(t,r)
        sigma_tr = g0_tr.copy()    # the 2 means second level of approximation, must be fermionic
        sigma_tr.data[..., 0, 0] = V_mtr.data[..., 0, 0, 0, 0] * g0_tr.data[..., 0, 0]
        
        # transform Sigma(t,r) to Sigma(w,k)
        #self.Sigma2_dlr_wk = fourier_tr_to_wk(Sigma2_dlr_tr)
        sigma_wk = fourier_tr_to_wk(sigma_tr)

        return sigma_wk
        
    def _calc_sigma_TPSC_dynamic(self):
        """
        Calculates the dynamic part of the second-level TPSC approximation of the self-energy.
        """

        # define effective potential
        V_wk = -self.U/8*(3*self.Usp*self.chisp_wk + self.Uch*self.chich_wk)

        V_wr = fourier_wk_to_wr(V_wk)
        V_tr = fourier_wr_to_tr(V_wr)

        g0_wr = fourier_wk_to_wr(self.g0_wk)
        g0_tr = fourier_wr_to_tr(g0_wr)
        
        sigma_tr = gw_dynamic_sigma(V_tr, g0_tr)

        sigma_wr = fourier_tr_to_wr(sigma_tr)
        sigma_wk = fourier_wr_to_wk(sigma_wr)

        if False:
            sigma_wk_ref = self._calc_sigma_deprecated()
            print('--> sigma test')
            np.testing.assert_array_almost_equal(sigma_wk.data, sigma_wk_ref.data)
            print('ok!')

        return sigma_wk
        
    def _get_sigma(self):
        """
        Calculates the self-energy including the Hartree-Fock term.
        """

        sigma_wk_TPSC_dynamic = self._calc_sigma_TPSC_dynamic()
        self.sigma_wk = sigma_wk_TPSC_dynamic.copy()
        self.sigma_wk.data[:] += self.U*self.n/2

### do not use, under construction ###
    def _imtime_bubble_chi2_wk(self):
        """
        Calculates chi2(r,tau) = - G2(r,tau)*G0(-r,-tau) - G2(-r,-tau)*G0(r,tau) in wk-space.

        Requires
        --------
        self.g0_dlr_wk and self.g2_dlr_wk must have been calculated

        Parameters
        ----------
        self                :   self

        Returns
        -------
        self.chi2_dlr_wk    :   second-level approximation of bubble in TPSC
        """

        # Fourier transform Gs to real space
        g2_tr = fourier_wk_to_tr(self.g_wk)
        g2_mtr = fourier_wk_to_mtr(self.g_wk)
        g0_tr = fourier_wk_to_tr(self.g0_wk)
        g0_mtr = fourier_wk_to_mtr(self.g0_wk)

        # calculate chi2
        chi2_tr = fourier_wk_to_tr(self.chi0_wk).copy()
        chi2_tr.data[:,0,0,0,0] = -g2_tr.data[:,0,0]*g0_mtr.data[:,0,0] - g2_mtr.data[:,0,0]*g0_tr.data[:,0,0]

        self.chi2_wk = fourier_tr_to_wk(chi2_tr)

        if False:
            from triqs_tprf.lattice import chi0_tr_from_grt_PH

            chi0_1 = chi0_tr_from_grt_PH(g2_tr, g0_tr)
            chi0_2 = chi0_tr_from_grt_PH(g0_tr, g2_tr)

            chi2_tr_ref = chi0_1 + chi0_2

            np.testing.assert_array_almost_equal(chi2_tr_ref.data, chi2_tr.data)


    def __skip_keys(self):
        return []


    def __eq__(self, obj):

        if obj.__dict__.keys() != self.__dict__.keys():
            return False

        for key in self.__dict__.keys():
            if key not in self.__skip_keys():
                a = getattr(self, key)
                b = getattr(obj, key)
                if not np.equal(a, b).all():
                    if type(a) == Gf and np.isclose(a.data, b.data).all(): continue
                    return False

        return True


    def __reduce_to_dict__(self):
        d = self.__dict__.copy()
        keys = set(d.keys()).intersection(self.__skip_keys())
        for key in keys: del d[key]
        return d


    @classmethod
    def __factory_from_dict__(cls, name, d):
        arg_keys = ['n', 'U', 'wmesh', 'e_k']
        argv_keys = ['docc', 'verbose']
        verbose = d['verbose']
        d['verbose'] = False # -- Suppress printouts on reconstruction from dict
        ret = cls(*[ d[key] for key in arg_keys ],
                  **{ key : d[key] for key in argv_keys })
        ret.__dict__.update(d)
        ret.verbose = verbose
        return ret


# -- Register Solver in Triqs formats

from h5.formats import register_class
register_class(tpsc_solver)

