/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include "../types.hpp"

namespace tprf {
  
/** Construct a non-interacting Matsubara frequency lattice Green's function $G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n)$

  Computes

  .. math::
     G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) = \left[ 
         (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
	 \right]^{-1}_{a\bar{b}},

  using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, chemical potential $\mu$,
  and a Matsubara frequency Green's function mesh.

  @param mu chemical potential :math:`\mu`
  @param ek discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
  @param mesh imaginary frequency mesh
  @return Matsubara frequency lattice Green's function $G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n)$
*/
gk_iw_t lattice_dyson_g0_wk(double mu, ek_vt e_k, g_iw_t::mesh_t mesh);
  
/** Construct an interacting Matsubara frequency lattice Green's function $G_{a\bar{b}}(\mathbf{k}, i\omega_n)$
   
 Computes

 .. math::
    G_{a\bar{b}}(\mathbf{k}, i\omega_n) = \left[ 
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a momentum independent Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n)$.

 @param mu chemical potential :math:`\mu`
 @param ek discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(\mathbf{k}, i\omega_n)$
 */
gk_iw_t lattice_dyson_g_wk(double mu, ek_vt e_k, g_iw_vt sigma_w);
g_iw_t lattice_dyson_g_w(double mu, ek_vt e_k, g_iw_vt sigma_w);

/** Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

    Computes: $G_{a\bar{b}}(\mathbf{r}, i\omega_n) = \mathcal{F}^{-1} \left\{ G_{a\bar{b}}(\mathbf{k}, i\omega_n) \right\}$

    @param gk k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\omega_n)`
    @return real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(\mathbf{r}, i\omega_n)`
 */
gr_iw_t fourier_wk_to_wr(gk_iw_vt g_wk);
  
/** Fast fourier transform of imaginary frequency Green's function from real-space to k-space

    Computes: $G_{a\bar{b}}(\mathbf{k}, i\omega_n) = \mathcal{F} \left\{ G_{a\bar{b}}(\mathbf{r}, i\omega_n) \right\}$

    @param gr k--space imaginary frequency Green's function :math:`G_{a\bar{b}}(\mathbf{r}, i\omega_n)`
    @return k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\omega_n)`
 */
gk_iw_t fourier_wr_to_wk(gr_iw_vt g_wr);

gr_tau_t fourier_wr_to_tr(gr_iw_vt g_wr, int nt=-1);

} // namespace tprf
