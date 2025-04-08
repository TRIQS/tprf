# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2025 by Xaver Landerl, Hugo U.R. Strand
# Authors: Xaver Landerl, Hugo U.R. Strand
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


class tpsc_solver:
    r"""
    Two-particle self-consistency [1] solver for single-band Hubbard models

    TPSC assumes RPA-like charge and spin susceptibilities
    .. math::
        \chi_{ch}(k) = \frac{\chi_0(k)}{1 + \frac{U_{ch}}{2}\chi_0(k)}, \quad
        \chi_{sp}(k) = \frac{\chi_0(k)}{1 - \frac{U_{sp}}{2}\chi_0(k)}
    with renormalized vertices and determines them such that the sum rules
    .. math::
        \frac{T}{N}\sum_{k}{\chi_{ch}(k)} = n + 2\braket{n_\uparrow n_\downarrow} - n^2, \quad
        \frac{T}{N}\sum_{k}{\chi_{sp}(k)} = n - 2\braket{n_\uparrow n_\downarrow}
    are fulfilled. This is done by imposing the Ansatz
    .. math::
        U_{sp}\braket{n_\uparrow}\braket{n_\downarrow} = U\braket{n_\uparrow n_\downarrow}.

    An approximation to the self-energy is obtained through
    .. math::
        \Sigma_\sigma(k) = Un_{-\sigma} + \frac{U}{8}\frac{T}{N}\sum_{q}{
        \left[3U_{sp}\chi_{sp}(k) + U_{ch}\chi_{ch}(k)\right]G_{0\sigma}(k+q)}.
        
    Determining the single-particle Green's function and self-energy,
    as well as the charge and spin susceptibilities.

    [1] Y.M. Vilk, A.-M.S. Tremblay, J. Phys. I France, v7, no 11, pp 1309-1368 (1997)
    https://doi.org/10.1051/jp1:1997135 and https://arxiv.org/abs/cond-mat/9702188
    """

    def __init__(self, n, U, wmesh, e_k, docc='tpsc_ansatz', verbose=True):
        r"""
        Initialize the TPSC solver.

        Parameters
        ----------
        n : double
            total electron density (spin up and down combined)
        U : double
            Hubbard interaction
        wmesh : MeshImFreq or MeshDLRImFreq
                imaginary frequency mesh
        e_k : Gf
              dispersion relation (excluding spin DOFs)
        docc : double, optional
               double occupancy (default: `tpsc_ansatz` using the TPSC-Ansatz [1])
        verbose : bool, optional
                  Verbose printouts (default: `True`)
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
        if self.use_tpsc_ansatz:
            self.vprint("Using the TPSC-Ansatz for the double occupancy.")
        else:
            self.vprint(f'docc = {self.docc} (not using the TPSC-ansatz)')
        self.vprint()
        self.vprint(f'len(wmesh) = {len(wmesh)}')
        self.vprint(f'nk = {len(self.e_k.mesh)}')
        self.vprint(f'beta = {self.beta}')
        self.vprint(f'n = {self.n}')
        self.vprint(f'U = {self.U}')
        self.vprint()


    def solve(self, calc_sigma=False, calc_g=False,
              Usp_tol=2e-12, Uch_tol=2e-12, Uch_max=100., Usp_epsilon=1e-7, chi0_wk=None):
        r"""
        Run the TPSC-calculation on the given model determining the screened
        spin and charge vertices `Usp` and `Uch`, respectively.

        Parameters
        ----------
        calc_sigma : bool, optional
                     Enable the TPSC-self-energy calculation (default: `False`)
        calc_g : bool, optional
                 Enable the lattice Green's function calculation (default: `False`)
        
        Usp_tol : double, optional
                  Tolerance for spin vertex `Usp` solution (default: 2e-12)
        Uch_tol : double, optional
                  Tolerance for charge vertex `Uch` solution (default: 2e-12)
        Uch_max : double, optional
                  maximum value of the charge vertex `Uch` to consider
                  in numerical search (default: 100.)
        Usp_epsilon : double, optional
                Offset from maximum value of the spin vertex `Usp` to consider
                Increases numerical stability of root search (default: 1e-7)
        chi0_wk : Gf, optional
                  Enables passing of a partially dressed susceptibility (default: None)
        """

        if calc_g == True: calc_sigma = True

        if chi0_wk is not None:
            self.chi0_wk = chi0_wk
        else:
            self.vprint("--> lattice_dyson_g0_wk")
            self.g0_wk = self._calc_g0_wk(self.n)
            nw = -1 if isinstance(self.wmesh, MeshDLRImFreq) else self.wmesh.n_iw # tprf api FIXME!

            self.vprint("--> imtime_bubble_chi0_wk")
            self.chi0_wk = 2*imtime_bubble_chi0_wk(self.g0_wk, nw=nw, verbose=False)

        self.vprint("--> get Usp")
        self.Usp = self._solve_Usp(self.chi0_wk, Usp_tol, Usp_epsilon)
        
        self.vprint("--> get Uch")
        self.Uch = self._solve_Uch(self.chi0_wk, Uch_tol, Uch_max)

        self.vprint("--> get chisp_wk and chich_wk")
        self.chisp_wk = self._solve_rpa(self.chi0_wk, +0.5 * self.Usp)
        self.chich_wk = self._solve_rpa(self.chi0_wk, -0.5 * self.Uch)

        chi_sum_rule = np.abs(self._get_density(self.chisp_wk + self.chich_wk) - (2*self.n - self.n**2))
        assert( chi_sum_rule < np.min([Usp_tol, Uch_tol]) )
        
        if self.use_tpsc_ansatz == True:
            self.vprint("--> get double occupancy")
            self.docc = self.Usp/self.U*self.n*self.n/4
        
        # print out results
        self.vprint()
        self.vprint("Summary first level approximation:")
        self.vprint(f"    Usp = {self.Usp}, Uch = {self.Uch}")
        self.vprint(f"    <n_up*n_down> = {self.docc}")
        self.vprint()

        if calc_sigma:
            self.vprint("--> get sigma_wk")
            self.sigma_wk = self.get_sigma()
        
        if calc_g:
            self.vprint("--> get g_wk")
            self.g_wk, self.mu = self._calc_g_wk(self.n, self.sigma_wk)

            self.vprint()
            self.vprint('Check sum rules of second-level approximation:')

            tr_SG0 = self._get_density(self.sigma_wk * self.g0_wk)
            tr_SG  = self._get_density(self.sigma_wk * self.g_wk)
            rel_diff = np.abs((tr_SG0 - tr_SG) / tr_SG0)

            self.vprint(f'    Tr[Sigma * G0] = {tr_SG0 - self.U * self.docc}')
            self.vprint(f'    Tr[Sigma * G ] = {tr_SG  - self.U * self.docc}')
            self.vprint(f'    Relative difference = {rel_diff:2.2E}')                    


    def _solve_Usp(self, chi0_wk, Usp_tol, Usp_epsilon):
        r"""
        Calculates screened spin vertex Usp given a bare susceptibility.

        Parameters
        ----------
        chi0_wk : Gf
                  bare susceptibility
        Usp_tol : double
                  Tolerance for spin vertex `Usp` solution
        Usp_epsilon : double
                      Offset from maximum value of the spin vertex `Usp` to consider
                      Increases numerical stability of root search

        Returns
        -------
        Usp : double
              Screened spin vertex `Usp`
        """

        def Usp_root(Usp):
            chi_wk = self._solve_rpa(chi0_wk, 0.5 * Usp)
            tr_chi = self._get_density(chi_wk)

            if self.use_tpsc_ansatz == True:
                RHS = self.n - Usp/self.U*self.n*self.n/2
            else:
                RHS = self.n - 2*self.docc

            return tr_chi - RHS
        
        # here chi_sp diverges (Usp_epsilon for numerical stability)
        Usp_max = 2.0/np.amax(self.chi0_wk.data).real - Usp_epsilon

        Usp = brentq(Usp_root, 0.0, Usp_max, xtol=Usp_tol)
        return Usp


    def _solve_Uch(self, chi0_wk, Uch_tol, Uch_max):
        r"""
        Calculates screened charge vertex `Uch` given a bare susceptibility.

        Parameters
        ----------
        chi0_wk : Gf
                  bare susceptibility
        Uch_tol : double
                  Tolerance for spin vertex `Uch` solution
    	Uch_max : double
                  maximum value of the charge vertex `Uch` to consider
                  in numerical search

        Returns
        -------
        Uch : double
              Screened spin vertex `Uch`
        """

        def Uch_root(Uch):
            chi_wk = self._solve_rpa(chi0_wk, -0.5 * Uch)
            tr_chi = self._get_density(chi_wk)
            
            if self.use_tpsc_ansatz == True:
                RHS = self.n + self.Usp/self.U*self.n*self.n/2 - self.n*self.n
            else:
                RHS = self.n + 2*self.docc - self.n*self.n
            
            return tr_chi - RHS
        
        Uch = brentq(Uch_root, 0.0, Uch_max, xtol=Uch_tol)
        return Uch


    def _solve_rpa(self, chi0_wk, U_vert):
        r"""
        Calculates an rpa-susceptibility with a scalar vertex from the non-interacting susceptibility.

        Parameters
        U_vert : double
                 RPA vertex

        Returns
        -------
        chi_wk : Gf
                 RPA-susceptibility with given vertex
        """

        if False:
            # TPRF RPA routine (its tensor/matrix products gives a large overhead for the scalar case)
            V = U_vert * np.ones((1, 1, 1, 1), dtype=complex)
            chi_wk = solve_rpa_PH(chi0_wk, V)

        chi_wk = chi0_wk.copy()
        chi_wk.data[:] = chi0_wk.data/(1 - U_vert * chi0_wk.data) # FIXME: 1/2 fator to U_vert
        return chi_wk

    
    def _get_density(self, g_wk):
        r"""
        Calculates
        .. math::
            \sum_{k}{g(k)}

        Parameters
        ----------
        g_wk : Gf
               Green's function in Matsubara frequency and momentum space
            
        Returns
        -------
        dens : double
               Density of passed Green's function
        """

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
        r"""
        Calculates the non-interacting Green's function
        .. math::
            g_0(k) = (i\omega_n + \mu - \varepsilon(\mathbf{k}))^{-1}
        using the solvers dispersion relation.
        Chemical potential is chosen such that g0_wk has the correct density.

        Parameters
        ----------
        target_density : double
                         total electron density (spin up and down combined)
        
        Returns
        -------
        g0_wk : Gf
                non-interacting Green's function
        """

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
        r"""
        Calculates the interacting Green's function
        .. math::
            g_\sigma(k) = (i\omega_n + \mu - \varepsilon(\mathbf{k}) - \Sigma_\sigma(k))
        using the solvers dispersion relation.
        Chemical potential is chosen such that g0_wk has the correct density.

        Parameters
        ----------
        target_density : double
                         total electron density (spin up and down combined)
        sigma_wk : Gf
                   self-energy
        
        Returns
        -------
        g_wk : Gf
               non-interacting Green's function
        mu : double
             chemical potential
        """

        # find the mu that leads to the correct density
        mu_min, mu_max = np.min(self.e_k.data.real), np.max(self.e_k.data.real)

        def target_function(mu):
            g_wk = lattice_dyson_g_wk(mu, self.e_k, sigma_wk)
            density = self._get_density(g_wk)
            return target_density/2 - density
    
        mu = brentq(target_function, mu_min, mu_max)
        g_wk = lattice_dyson_g_wk(mu, self.e_k, sigma_wk)
        
        return g_wk, mu
    
    
    def get_sigma_tpsc_dynamic_numpy(self):
        r"""
        Calculates the dynamic part of the second-level TPSC approximation of the self-energy.
        .. math::
            \Sigma_\sigma^{dyn, TPSC}(k) = \frac{U}{8}\frac{T}{N}\sum_{q}{
            \left[3U_{sp}\chi_{sp}(k) + U_{ch}\chi_{ch}(k)\right]G_{0\sigma}(k+q)}

        Uses numpy-arrays instead of built-in TRIQS functionality.

        Returns
        -------
        sigma_wk : Gf
                   dynamic TPSC-self energy
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
        

    def get_sigma_tpsc_dynamic(self):
        r"""
        Calculates the dynamic part of the second-level TPSC approximation of the self-energy.
        .. math::
            \Sigma_\sigma^{dyn, TPSC}(k) = \frac{U}{8}\frac{T}{N}\sum_{q}{
            \left[3U_{sp}\chi_{sp}(k) + U_{ch}\chi_{ch}(k)\right]G_{0\sigma}(k+q)}

        Returns
        -------
        sigma_wk : Gf
                   dynamic TPSC-self energy
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


    def get_sigma(self):
        r"""
        Calculates the full second-level TPSC approximation to the self-energy.
        .. math::
            \Sigma_\sigma^{TPSC}(k) = Un_{-\sigma} \frac{U}{8}\frac{T}{N}\sum_{q}{
            \left[3U_{sp}\chi_{sp}(k) + U_{ch}\chi_{ch}(k)\right]G_{0\sigma}(k+q)}
        
        Returns
        -------
        sigma_wk : Gf
                   TPSC-self energy
        """

        sigma_dyn_wk = self.get_sigma_tpsc_dynamic()
        sigma_wk = sigma_dyn_wk.copy()
        sigma_wk.data[:] += self.U*self.n/2
        return sigma_wk


    def get_GG0_bubble_chi2_wk(self):
        r"""
        Calculates the partially dressed susceptibility
        .. math::
            chi_2(k) = -\frac{T}{N}\sum_{q}{\left[
            G_\sigma^{TPSC}(q)G_{0\sigma}(q+k) + G_\sigma^{TPSC}(q+k)G_{0\sigma}(q)\right]}
        used in the TPSC+ calculation [2]

        [2] C. Gauvin-Ndiaye, C. Lahaie, Y.M. Vilk, A.-M.S. Tremblay Phys.Rev.B 108, 075144, 2023
        https://doi.org/10.1103/PhysRevB.108.075144 and https://doi.org/10.48550/arXiv.2305.19219
        
        Returns
        -------
        chi2_wk : Gf
        """

        # Fourier transform Gs to real space
        g2_tr = fourier_wk_to_tr(self.g_wk)
        g0_tr = fourier_wk_to_tr(self.g0_wk)

        from triqs_tprf.lattice import chi0_tr_from_grt_PH

        chi2_tr = \
            chi0_tr_from_grt_PH(g2_tr, g0_tr) + \
            chi0_tr_from_grt_PH(g0_tr, g2_tr)

        chi2_wk = fourier_tr_to_wk(chi2_tr)

        if False:
            g2_mtr = fourier_wk_to_mtr(self.g_wk)
            g0_mtr = fourier_wk_to_mtr(self.g0_wk)

            # calculate chi2
            chi2_tr_ref = fourier_wk_to_tr(self.chi0_wk).copy()
            chi2_tr_ref.data[:, :, 0, 0, 0, 0] = \
                -  g2_tr.data[:, :, 0, 0] * g0_mtr.data[:, :, 0, 0] \
                - g2_mtr.data[:, :, 0, 0] *  g0_tr.data[:, :, 0, 0]

            np.testing.assert_array_almost_equal(chi2_tr_ref.data, chi2_tr.data)

        return chi2_wk


    def vprint(self, string=None):
        if self.verbose == True:
            if string is not None:
                print(string)
            else:
                print()


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

    Parameters
    ----------
    g_wk : Gf
           Green's function in wk representation

    Returns
    -------
    g_tr : Gf
           Fourier-transform of passed Green's function
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
    Fourier-transforms a Green's function from wk to mtr representation.

    Parameters
    ----------
    g_wk : Gf
           Green's function in wk representation

    Returns
    -------
    g_mtr : Gf
            Fourier-transform of passed Green's function
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
    Fourier-transforms a Green's function from tr to wk representation.

    Parameters
    ----------
    g_tr : Gf
           Green's function in wk representation

    Returns
    -------
    g_wk : Gf
           Fourier-transform of passed Green's function
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
