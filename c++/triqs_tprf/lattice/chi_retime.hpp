/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, H. U.R. Strand
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

std::tuple<g_Tk_t, g_Tk_t> g0_Tk_les_gtr_from_e_k(e_k_cvt e_k, gf_mesh<retime> Tmesh, double beta);
  
/** Generalized susceptibility real time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r})`

  Computes

  .. math::
     \chi^{(0)R}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
     i G^<_{d\bar{a}}(t, \mathbf{r}) G^>_{b\bar{c}}(-t, -\mathbf{r})
     - i G^>_{d\bar{a}}(t, \mathbf{r}) G^<_{b\bar{c}}(-t, -\mathbf{r})

  @param g_Tr_les Lesser real time Green's function in real-space, :math:`G^<_{a\bar{b}}(t, \mathbf{r})`.
  @param g_Tr_gtr Greater real time Green's function in real-space, :math:`G^>_{a\bar{b}}(t, \mathbf{r})`.
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r})` in real time and real-space.
 */
chi_Tr_t chi0_Tr_from_g_Tr_PH(g_Tr_cvt g_Tr_les, g_Tr_cvt g_Tr_gtr);

} // namespace triqs_tprf
