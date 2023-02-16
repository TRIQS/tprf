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

namespace triqs_tprf {

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
  g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, mesh::imfreq mesh);

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
  g_Dwk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, mesh::dlr_imfreq mesh);

  /** Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
   
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

  /** Construct a non-interacting real frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(\omega, \mathbf{k})`

  Computes

  .. math::
     G^{(0)}_{a\bar{b}}(\omega, \mathbf{k}) = \left[
     (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
     \right]^{-1}_{a\bar{b}},

  using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, chemical potential $\mu$,
  broadening $\delta$, and a real frequency Green's function mesh.

  @param mu chemical potential :math:`\mu`
  @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
  @param mesh real frequency mesh
  @param delta broadening :math:`\delta`
  @return Matsubara frequency lattice Green's function $G^{(0)}_{a\bar{b}}(\omega, \mathbf{k})$
*/
  g_fk_t lattice_dyson_g0_fk(double mu, e_k_cvt e_k, mesh::refreq mesh, double delta);

  /** Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_wk imaginary frequency self-energy $\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})$
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
  g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_wk_cvt sigma_wk);

  /** Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_wk imaginary frequency self-energy $\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})$
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
  g_Dwk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_Dwk_cvt sigma_wk);

  /** Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n)$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_w imaginary frequency self-energy $\Sigma_{\bar{a}b}(i\omega_n)$
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
  g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w);

  /** Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
   
 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a Matsubara frequency 
 self energy $\Sigma_{\bar{a}b}(i\omega_n)$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_w imaginary frequency self-energy $\Sigma_{\bar{a}b}(i\omega_n)$
 @return Matsubara frequency lattice Green's function $G_{a\bar{b}}(i\omega_n, \mathbf{k})$
 */
  g_Dwk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_Dw_cvt sigma_w);

  /** Construct an interacting real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
   
 Computes

 .. math::
    G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega, \mathbf{k})
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, broadening $\delta$, and a real frequency 
 self energy $\Sigma_{\bar{a}b}(\omega, \mathbf{k})$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_fk real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega, \mathbf{k})`
 @param delta broadening :math:`\delta`
 @return real frequency lattice Green's function $G_{a\bar{b}}(\omega, \mathbf{k})$
 */
  g_fk_t lattice_dyson_g_fk(double mu, e_k_cvt e_k, g_fk_cvt sigma_fk, double delta);

  /** Construct an interacting real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
   
 Computes

 .. math::
    G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, broadening $\delta$, and a real frequency 
 self energy $\Sigma_{\bar{a}b}(\omega)$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_fk real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega)`
 @param delta broadening :math:`\delta`
 @return real frequency lattice Green's function $G_{a\bar{b}}(\omega, \mathbf{k})$
 */
  g_fk_t lattice_dyson_g_fk(double mu, e_k_cvt e_k, g_f_cvt sigma_f, double delta);

  /** Construct an interacting Matsubara frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function :math:`G_{a\bar{b}}(i\omega_n)`
   
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
  g_Dw_t lattice_dyson_g_w(double mu, e_k_cvt e_k, g_Dw_cvt sigma_w);

  /** Construct an interacting real frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function :math:`G_{a\bar{b}}(\omega)`
   
 Computes

 .. math::
    G_{a\bar{b}}(\omega) = \frac{1}{N_k} \sum_\mathbf{k} \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
 chemical potential $\mu$, and a momentum independent real frequency 
 self energy $\Sigma_{\bar{a}b}(\omega)$.

 @param mu chemical potential :math:`\mu`
 @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
 @param sigma_f real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega)`
 @param delta broadening :math:`\delta`
 @return Real frequency lattice Green's function $G_{a\bar{b}}(\omega, \mathbf{k})$
 */
  g_f_t lattice_dyson_g_f(double mu, e_k_cvt e_k, g_f_cvt sigma_f, double delta);

  /** Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(i\omega_n, \mathbf{k})\right\}`

    @param g_wk k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
    @return real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
 */
  g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk);
  g_Dwr_t fourier_wk_to_wr(g_Dwk_cvt g_wk);

  /** Fast fourier transform of imaginary frequency Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

    @param g_wr real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
    @return k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
 */
  g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr);
  g_Dwk_t fourier_wr_to_wk(g_Dwr_cvt g_wr);

  /** Fast fourier transform of real-space Green's function from Matsubara frequency to imaginary time

    Computes: :math:`G_{a\bar{b}}(\tau, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

    @param g_wr real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
    @return real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`
 */
  g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt = -1);
  g_Dtr_t fourier_wr_to_tr(g_Dwr_cvt g_wr, int nt = -1);

  /** Fast fourier transform of real-space Green's function from imaginary time to Matsubara frequency

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(\tau, \mathbf{r}) \right\}`

    @param g_tr real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`
    @return real-space Matsubara frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
 */
  g_wr_t fourier_tr_to_wr(g_tr_cvt g_tr, int nw = -1);
  g_Dwr_t fourier_tr_to_wr(g_Dtr_cvt g_tr, int nw = -1);
  
  /** Inverse fast fourier transform of real frequency Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(\omega, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(\omega, \mathbf{k})\right\}`

    @param g_fk k-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
    @return real-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{r})`
 */
  g_fr_t fourier_fk_to_fr(g_fk_cvt g_fk);

  /** Fast fourier transform of real frequency Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(\omega, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(\omega, \mathbf{r}) \right\}`

    @param g_fr real-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{r})`
    @return k-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
 */
  g_fk_t fourier_fr_to_fk(g_fr_cvt g_fr);

  /** Inverse fast fourier transform of real time Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(t, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(t, \mathbf{k})\right\}`

    @param g_Tk k-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{k})`
    @return real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`
 */
  g_Tr_t fourier_Tk_to_Tr(g_Tk_cvt g_Tk);

  /** Fast fourier transform of real time Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(t, \mathbf{r}) \right\}`

    @param g_Tr real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`
    @return k-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{k})`
 */
  g_Tk_t fourier_Tr_to_Tk(g_Tr_cvt g_Tr);

  /** Fast fourier transform of real time Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}c\bar{d}}(t, \mathbf{r}) \right\}`

    @param g_Tr real-space real time Green's function :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{r})`
    @return k-space real time Green's function :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{k})`
 */
  chi_Tk_t fourier_Tr_to_Tk(chi_Tr_cvt chi_Tr);

  /** Fast fourier transform of real time Green's function from k-space to real-space

    Computes: :math:`G_{a\bar{b}}(t, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(t, \mathbf{k})\right\}`

    @param g_Tk k-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{k})`
    @return real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`
 */
  chi_Tr_t fourier_Tk_to_Tr(chi_Tk_cvt chi_Tk);
  
} // namespace triqs_tprf
