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
  
gk_iw_t g0k_from_ek(double mu, ek_cvt ek, g_iw_t::mesh_t mesh);
gk_iw_t gk_from_ek_sigma(double mu, ek_cvt ek, g_iw_cvt sigma);

gr_iw_t gr_from_gk(gk_iw_cvt gk);
gk_iw_t gk_from_gr(gr_iw_cvt gr, brillouin_zone const &bz);

chi0r_t chi0r_from_gr_PH(int nw, int nnu, gr_iw_cvt gr);

chi0r_t chi0r_from_chi0q(chi0q_cvt chi0q);
chi0q_t chi0q_from_chi0r(chi0r_cvt chi0r, brillouin_zone const & bz);

} // namespace tprf
