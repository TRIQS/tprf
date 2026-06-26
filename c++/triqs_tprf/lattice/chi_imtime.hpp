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

/** Generalized susceptibility imaginary time bubble in the particle-hole channel \f$ \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$

  Computes

  \f[
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
     - G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})
  \f]

  @param g_tr Imaginary time Green's function in real-space, \f$ G_{a\bar{b}}(\tau, \mathbf{r}) \f$.
  @return Generalized susceptibility \f$ \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$ in imaginary time and real-space.
 */
chi_tr_t chi0_tr_from_grt_PH(g_tr_cvt g_tr);
chi_tr_t chi0_tr_from_grt_PH(g_tr_cvt g_tr, g_tr_cvt g_bwd_tr);
chi_Dtr_t chi0_tr_from_grt_PH(g_Dtr_cvt g_tr);
chi_Dtr_t chi0_tr_from_grt_PH(g_Dtr_cvt g_tr, g_Dtr_cvt g_bwd_tr);
chi_wr_t chi0_wr_from_grt_PH(g_tr_cvt g_tr, int nw);
chi_wr_t chi0_wr_from_grt_PH(g_tr_cvt g_tr, g_tr_cvt g_bwd_tr, int nw);

/** Generalized susceptibility zero imaginary frequency bubble in the particle-hole channel \f$ \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) \f$

  Computes

  \f[
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r}) =
     - \int_0^\beta d\tau \,
     G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})
  \f]

  @param g_tr Imaginary time Green's function in real-space, \f$ G_{a\bar{b}}(\tau, \mathbf{r}) \f$.
  @return Generalized susceptibility \f$ \chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r}) \f$ in real-space.
 */
chi_wr_t chi0_w0r_from_grt_PH(g_tr_cvt g_tr);
chi_wr_t chi0_w0r_from_grt_PH(g_tr_cvt g_tr, g_tr_cvt g_bwd_tr);
chi_wr_t chi0_w0r_from_grt_PH(g_Dtr_cvt g_tr);
chi_wr_t chi0_w0r_from_grt_PH(g_Dtr_cvt g_tr, g_Dtr_cvt g_bwd_tr);

/** Static susceptibility calculation \f$ \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) \f$
   
  Explicit calculation of the static, zero frequency response, by 2nd order trapetzoidal 
  integration in imaginary time, i.e.
  
  \f[
     \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) =
         \int_0^\beta d\tau \, \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) 
  \f]

  @param chi_tr Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$ 
                in imaginary time and real space.
  @return Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) \f$ 
          at zero Matsubara frequency and real-space.
 */
chi_wr_t chi_w0r_from_chi_tr(chi_tr_cvt chi_tr);
  
/** Parallel Fourier transform from  \f$ \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$ to \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$

  Computes

  \f[
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \mathcal{F}_{\tau \rightarrow \omega} \left\{
         \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) 
         \right\}
  \f]

  @param chi_tr Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$ 
                in imaginary time and real space.
  @return Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$ 
          in Matsubara frequency and real-space.
 */
chi_wr_t chi_wr_from_chi_tr(chi_tr_cvt chi_tr, int nw);
chi_Dwr_t chi_wr_from_chi_tr(chi_Dtr_cvt chi_tr, int nw);

/** Fourier transform from \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$ to \f$ \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$

  Computes

  \f[
         \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
         \mathcal{F}_{\omega \rightarrow \tau} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \right\}
  \f]

  @param chi_tr Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) \f$ 
                in imaginary time and real space.
  @return Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$ 
          in Matsubara frequency and real-space.
 */
chi_tr_t chi_tr_from_chi_wr(chi_wr_cvt chi_wr, int ntau=-1);
chi_Dtr_t chi_tr_from_chi_wr(chi_Dwr_cvt chi_wr, int ntau=-1);

/** Parallel Fourier transform from \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$ to \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) \f$

  Computes

  \f[
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
         \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{k}} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) 
         \right\}
  \f]

  @param chi_wr Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$ 
                in Matsubara frequency and real space.
  @return Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) \f$ 
          in Matsubara frequency and momentum space.
 */
chi_wk_t chi_wk_from_chi_wr(chi_wr_cvt chi_wr);
chi_Dwk_t chi_wk_from_chi_wr(chi_Dwr_cvt chi_wr);

/** Parallel Fourier transform from \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) \f$ to \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$

  Computes

  \f[
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \mathcal{F}_{\mathbf{k} \rightarrow \mathbf{r}} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) 
         \right\}
  \f]

  @param chi_wr Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) \f$ 
                in imaginary time and momentum space.
  @return Generalized susceptibility \f$ \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) \f$ 
          in Matsubara frequency and real space.
 */
chi_wr_t chi_wr_from_chi_wk(chi_wk_cvt chi_wk);
chi_Dwr_t chi_wr_from_chi_wk(chi_Dwk_cvt chi_wk);

target_value_t<chi_t_t>::regular_type chi_trapz_tau(chi_t_cvt chi_t);
target_value_t<chi_t_t>::regular_type integrate_dlr_tau(chi_Dt_cvt chi_t);

} // namespace triqs_tprf
