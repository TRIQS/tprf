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

namespace tprf {

  template <Channel_t CH> g2_iw_t chi0_from_gg2(g_iw_cvt g, g2_iw_cvt g2);
  template <Channel_t CH> g2_iw_t chi_from_gg2(g_iw_cvt g, g2_iw_cvt g2);

// ----------------------------------------------------

/*
void calculate_chi0_and_chi_ph(g2_iw_vt chi0, g2_iw_vt chi, g_iw_cvt g,
g2_iw_cvt g2) {
  double beta = g.mesh().domain().beta;

  chi0(Omega, n, np)(a, b, c, d) << - beta * kronecker(n, np) * g(n)(d, a) *
g(Omega + n)(b, c);

  chi(Omega, n, np)(a, b, c, d) <<
    g2(Omega, n, np)(a, b, c, d) - beta * kronecker(Omega) * g(n)(b, a) *
g(np)(d, c);
}
*/

/**

 */

/*
void calculate_chi0_and_chi_pp(g2_iw_vt chi0, g2_iw_vt chi, g_iw_cvt g,
g2_iw_cvt g2) {
  double beta = g.mesh().domain().beta;

  chi0(Omega, n, np)(a, b, c, d) << - beta * kronecker(n, np) * g(n)(d, a) *
g(Omega - n)(b, c);

  chi(Omega, n, np)(a, b, c, d) <<
    g2(Omega, n, np)(a, b, c, d) - beta * kronecker(n + np, Omega) * g(n)(b, a)
* g(np)(d, c);
}
*/

} // namespace tprf
