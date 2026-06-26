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

  /** Construct a non-interacting Matsubara frequency lattice Green's function \f$ G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$

  Computes

  \f[
     G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
     (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
     \right]^{-1}_{a\bar{b}},
  \f]

  using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, chemical potential \f$ \mu \f$,
  and a Matsubara frequency Green's function mesh.

  @param mu chemical potential \f$ \mu \f$
  @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
  @param mesh imaginary frequency mesh
  @return Matsubara frequency lattice Green's function \f$ G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
  */
  g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, mesh::imfreq mesh);

  /** Construct a non-interacting Matsubara frequency lattice Green's function \f$ G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$

  Computes

  \f[
     G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
     (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
     \right]^{-1}_{a\bar{b}},
  \f]

  using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, chemical potential \f$ \mu \f$,
  and a Matsubara frequency Green's function mesh.

  @param mu chemical potential \f$ \mu \f$
  @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
  @param mesh imaginary frequency mesh
  @return Matsubara frequency lattice Green's function \f$ G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
  */
  g_Dwk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, mesh::dlr_imfreq mesh);

  /** Construct an interacting Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a momentum independent Matsubara frequency 
 self energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_w imaginary frequency self-energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$
 @return Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w);

  /** Construct a non-interacting real frequency lattice Green's function \f$ G^{(0)}_{a\bar{b}}(\omega, \mathbf{k}) \f$

  Computes

  \f[
     G^{(0)}_{a\bar{b}}(\omega, \mathbf{k}) = \left[
     (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
     \right]^{-1}_{a\bar{b}},
  \f]

  using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, chemical potential \f$ \mu \f$,
  broadening \f$ \delta \f$, and a real frequency Green's function mesh.

  @param mu chemical potential \f$ \mu \f$
  @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
  @param mesh real frequency mesh
  @param delta broadening \f$ \delta \f$
  @return Matsubara frequency lattice Green's function \f$ G^{(0)}_{a\bar{b}}(\omega, \mathbf{k}) \f$
*/
  g_fk_t lattice_dyson_g0_fk(double mu, e_k_cvt e_k, mesh::refreq mesh, double delta);

  /** Construct an interacting Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a Matsubara frequency 
 self energy \f$ \Sigma_{\bar{a}b}(i\omega_n, \mathbf{k}) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_wk imaginary frequency self-energy \f$ \Sigma_{\bar{a}b}(i\omega_n, \mathbf{k}) \f$
 @return Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_wk_cvt sigma_wk);

  /** Construct an interacting Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a Matsubara frequency 
 self energy \f$ \Sigma_{\bar{a}b}(i\omega_n, \mathbf{k}) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_wk imaginary frequency self-energy \f$ \Sigma_{\bar{a}b}(i\omega_n, \mathbf{k}) \f$
 @return Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_Dwk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_Dwk_cvt sigma_wk);

  /** Construct an interacting Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a Matsubara frequency 
 self energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_w imaginary frequency self-energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$
 @return Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w);

  /** Construct an interacting Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a Matsubara frequency 
 self energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_w imaginary frequency self-energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$
 @return Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_Dwk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_Dw_cvt sigma_w);

  /** Construct an interacting real frequency lattice Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega, \mathbf{k})
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, broadening \f$ \delta \f$, and a real frequency 
 self energy \f$ \Sigma_{\bar{a}b}(\omega, \mathbf{k}) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_fk real frequency self-energy \f$ \Sigma_{\bar{a}b}(\omega, \mathbf{k}) \f$
 @param delta broadening \f$ \delta \f$
 @return real frequency lattice Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
 */
  g_fk_t lattice_dyson_g_fk(double mu, e_k_cvt e_k, g_fk_cvt sigma_fk, double delta);

  /** Construct an interacting real frequency lattice Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, broadening \f$ \delta \f$, and a real frequency 
 self energy \f$ \Sigma_{\bar{a}b}(\omega) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_fk real frequency self-energy \f$ \Sigma_{\bar{a}b}(\omega) \f$
 @param delta broadening \f$ \delta \f$
 @return real frequency lattice Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
 */
  g_fk_t lattice_dyson_g_fk(double mu, e_k_cvt e_k, g_f_cvt sigma_f, double delta);

  /** Construct an interacting Matsubara frequency local (\f$ \mathbf{r}=\mathbf{0} \f$) lattice Green's function \f$ G_{a\bar{b}}(i\omega_n) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(i\omega_n) = \frac{1}{N_k} \sum_\mathbf{k} \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a momentum independent Matsubara frequency 
 self energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_w imaginary frequency self-energy \f$ \Sigma_{\bar{a}b}(i\omega_n) \f$
 @return Matsubara frequency lattice Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_w_t lattice_dyson_g_w(double mu, e_k_cvt e_k, g_w_cvt sigma_w);
  g_Dw_t lattice_dyson_g_w(double mu, e_k_cvt e_k, g_Dw_cvt sigma_w);

  /** Construct an interacting real frequency local (\f$ \mathbf{r}=\mathbf{0} \f$) lattice Green's function \f$ G_{a\bar{b}}(\omega) \f$
   
 Computes

 \f[
    G_{a\bar{b}}(\omega) = \frac{1}{N_k} \sum_\mathbf{k} \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
    \right]^{-1}_{a\bar{b}},
 \f]

 using a discretized dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, 
 chemical potential \f$ \mu \f$, and a momentum independent real frequency 
 self energy \f$ \Sigma_{\bar{a}b}(\omega) \f$.

 @param mu chemical potential \f$ \mu \f$
 @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
 @param sigma_f real frequency self-energy \f$ \Sigma_{\bar{a}b}(\omega) \f$
 @param delta broadening \f$ \delta \f$
 @return Real frequency lattice Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
 */
  g_f_t lattice_dyson_g_f(double mu, e_k_cvt e_k, g_f_cvt sigma_f, double delta);

  /** Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

    Computes: \f$ G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(i\omega_n, \mathbf{k})\right\} \f$

    @param g_wk k-space imaginary frequency Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
    @return real-space imaginary frequency Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \f$
 */
  g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk);
  g_Dwr_t fourier_wk_to_wr(g_Dwk_cvt g_wk);
  chi_Dwr_t fourier_wk_to_wr(chi_Dwk_cvt chi_wk);

  /** Fast fourier transform of imaginary frequency Green's function from real-space to k-space

    Computes: \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\} \f$

    @param g_wr real-space imaginary frequency Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \f$
    @return k-space imaginary frequency Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \f$
 */
  g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr);
  g_Dwk_t fourier_wr_to_wk(g_Dwr_cvt g_wr);
  chi_Dwk_t fourier_wr_to_wk(chi_Dwr_cvt chi_wr);

  /** Fast fourier transform of real-space Green's function from Matsubara frequency to imaginary time

    Computes: \f$ G_{a\bar{b}}(\tau, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\} \f$

    @param g_wr real-space imaginary frequency Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \f$
    @return real-space imaginary time Green's function \f$ G_{a\bar{b}}(\tau, \mathbf{r}) \f$
 */
  g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt = -1);
  g_Dtr_t fourier_wr_to_tr(g_Dwr_cvt g_wr, int nt = -1);
  chi_Dtr_t fourier_wr_to_tr(chi_Dwr_cvt chi_wr, int nt = -1);

  /** Fast fourier transform of real-space Green's function from imaginary time to Matsubara frequency

    Computes: \f$ G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(\tau, \mathbf{r}) \right\} \f$

    @param g_tr real-space imaginary time Green's function \f$ G_{a\bar{b}}(\tau, \mathbf{r}) \f$
    @return real-space Matsubara frequency Green's function \f$ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \f$
 */
  g_wr_t fourier_tr_to_wr(g_tr_cvt g_tr, int nw = -1);
  g_Dwr_t fourier_tr_to_wr(g_Dtr_cvt g_tr, int nw = -1);
  chi_Dwr_t fourier_tr_to_wr(chi_Dtr_cvt chi_tr, int nw = -1);
  
  /** Inverse fast fourier transform of real frequency Green's function from k-space to real space

    Computes: \f$ G_{a\bar{b}}(\omega, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(\omega, \mathbf{k})\right\} \f$

    @param g_fk k-space real frequency Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
    @return real-space real frequency Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{r}) \f$
 */
  g_fr_t fourier_fk_to_fr(g_fk_cvt g_fk);

  /** Fast fourier transform of real frequency Green's function from real-space to k-space

    Computes: \f$ G_{a\bar{b}}(\omega, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(\omega, \mathbf{r}) \right\} \f$

    @param g_fr real-space real frequency Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{r}) \f$
    @return k-space real frequency Green's function \f$ G_{a\bar{b}}(\omega, \mathbf{k}) \f$
 */
  g_fk_t fourier_fr_to_fk(g_fr_cvt g_fr);

  /** Inverse fast fourier transform of real time Green's function from k-space to real space

    Computes: \f$ G_{a\bar{b}}(t, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(t, \mathbf{k})\right\} \f$

    @param g_Tk k-space real time Green's function \f$ G_{a\bar{b}}(t, \mathbf{k}) \f$
    @return real-space real time Green's function \f$ G_{a\bar{b}}(t, \mathbf{r}) \f$
 */
  g_Tr_t fourier_Tk_to_Tr(g_Tk_cvt g_Tk);

  /** Fast fourier transform of real time Green's function from real-space to k-space

    Computes: \f$ G_{a\bar{b}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(t, \mathbf{r}) \right\} \f$

    @param g_Tr real-space real time Green's function \f$ G_{a\bar{b}}(t, \mathbf{r}) \f$
    @return k-space real time Green's function \f$ G_{a\bar{b}}(t, \mathbf{k}) \f$
 */
  g_Tk_t fourier_Tr_to_Tk(g_Tr_cvt g_Tr);

  /** Fast fourier transform of real time Green's function from real-space to k-space

    Computes: \f$ G_{a\bar{b}c\bar{d}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}c\bar{d}}(t, \mathbf{r}) \right\} \f$

    @param g_Tr real-space real time Green's function \f$ G_{a\bar{b}c\bar{d}}(t, \mathbf{r}) \f$
    @return k-space real time Green's function \f$ G_{a\bar{b}c\bar{d}}(t, \mathbf{k}) \f$
 */
  chi_Tk_t fourier_Tr_to_Tk(chi_Tr_cvt chi_Tr);
    
} // namespace triqs_tprf
