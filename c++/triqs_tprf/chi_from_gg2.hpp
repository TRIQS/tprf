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

#include "types.hpp"

namespace triqs_tprf {

template <Channel_t CH> g2_iw_t chi0_from_gg2(g_iw_cvt g, g2_iw_cvt g2);

/** Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Hole channel

     Computes

     .. math::
         \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega + \nu)

     @param g single particle Green's function :math:`G_{ab}(\nu)`
     @param g2 two-particle Green's function with the mesh to use for
   :math:`\chi^{(0)}`
     @return chi0 particle-hole bubble
   :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
  */
g2_iw_t chi0_from_gg2_PH(g_iw_vt g, g2_iw_vt g2);

/** Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Particle
   channel

     Computes

     .. math::
         \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega - \nu)

     @param g single particle Green's function :math:`G_{ab}(\nu)`
     @param g2 two-particle Green's function with the mesh to use for
   :math:`\chi^{(0)}`
     @return chi0 particle-particle bubble
   :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
  */
g2_iw_t chi0_from_gg2_PP(g_iw_vt g, g2_iw_vt g2);

template <Channel_t CH> g2_iw_t chi_from_gg2(g_iw_cvt g, g2_iw_cvt g2);

/** Generalized susceptibility :math:`\chi^{(0)} = G^{(2)} - GG` in the
   Particle-Hole channel

     Computes

     .. math::
         \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
         - \beta \delta_{\omega} G_{ba}(\nu) \cdot G_{dc}(\nu')

     @param g single particle Green's function :math:`G_{ab}(\nu)`
     @param g2 two-particle Green's function
   :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`
     @return chi generalized particle-hole susceptibility
   :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
  */
g2_iw_t chi_from_gg2_PH(g_iw_vt g, g2_iw_vt g2);

/** Generalized susceptibility :math:`\chi^{(0)} = G^{(2)} - GG` in the
   Particle-Particle channel

     Computes

     .. math::
         \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
         - \beta \delta_{\nu + \nu' - \omega} G_{ba}(\nu) \cdot G_{dc}(\nu')

     @param g single particle Green's function :math:`G_{ab}(\nu)`
     @param g2 two-particle Green's function
   :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`
     @return chi generalized particle-hole susceptibility
   :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
  */
g2_iw_t chi_from_gg2_PP(g_iw_vt g, g2_iw_vt g2);

/** Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Hole channel

     Computes

     .. math::
         \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau) =
         G_{da}(\nu) \cdot G_{bc}(\beta - \tay)

     @param g single particle Green's function :math:`G_{ab}(\tau)`
     @return chi0 particle-hole bubble
   :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau)`
  */
chi2_tau_t chi0_tau_from_g_tau_PH(g_tau_cvt g);
  
} // namespace triqs_tprf
