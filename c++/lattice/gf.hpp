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
  
/** Construct a non-interacting Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`

  Computes

  .. math::
     G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[ 
         (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
	 \right]^{-1}_{a\bar{b}},

  using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, chemical potential $\mu$,
  and a Matsubara frequency Green's function mesh.

  @param mu chemical potential :math:`\mu`
  @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
  @param mesh imaginary frequency mesh
  @return Matsubara frequency lattice Green's function $G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})$
*/
g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, gf_mesh<imfreq> mesh);
  
/** Construct an interacting Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[ 
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a momentum independent Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n)$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_w imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w);

/** Construct an interacting Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[ 
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a momentum independent Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_wk imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_wk_cvt sigma_wk);

/** Construct an interacting Matsubara frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function $G_{a\bar{b}}(i\omega_n)$
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n) = \frac{1}{N_k} \sum_\mathbf{k} \left[ 
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a momentum independent Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n)$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_w imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
g_w_t lattice_dyson_g_w(double mu, e_k_cvt e_k, g_w_cvt sigma_w);

/** Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \right\}`

    @param g_wk k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
    @return real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
 */
g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk);
  
/** Fast fourier transform of imaginary frequency Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

    @param g_wr real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
    @return k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
 */
g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr);

/** Fast fourier transform of real-space Green's function from Matsubara frequency to imaginary time

    Computes: :math:`G_{a\bar{b}}(\tau, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

    @param g_wr real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
    @return real-space imaginary time Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
 */
g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt=-1);

} // namespace tprf
