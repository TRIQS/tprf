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

/** Generalized susceptibility imaginary time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
     - G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})

  .. note::
     The imaginary time Green's function evaluator does not work together with
     triqs-clef expressions, see bug https://github.com/TRIQS/triqs/pull/558
     When the evaluator is fixed this routine should be reverted to the simpler
     clef expression (now commented out).

  @param grt Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real-space

 */
chi_tr_t chi0_tr_from_grt_PH(gr_tau_vt grt);

chi_wr_t chi0_w0r_from_grt_PH(gr_tau_vt grt);

/** Static susceptibility calculation :math:`\chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`
   
  Explicit calculation of the static, zero frequency response, by 2nd order trapetzoidal 
  integration in imaginary time, i.e.
  
  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) =
         \int_0^\beta d\tau \, \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) 

  @param chi_tr Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real space
  @return Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real-space
 */
chi_wr_t chi_w0r_from_chi_tr(chi_tr_vt chi_tr);
  
chi_wr_t chi_wr_from_chi_tr(chi_tr_vt chi_tr, int nw);

chi_wk_t chi_wk_from_chi_wr(chi_wr_vt chi_wr);

chi_t_t::zero_t chi_trapz_tau(chi_t_vt chi_t);

} // namespace tprf
