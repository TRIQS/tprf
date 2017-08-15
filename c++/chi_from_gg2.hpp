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

g2_iw_t chi0_from_gg2_PH(g_iw_vt g, g2_iw_vt g2);
g2_iw_t chi0_from_gg2_PP(g_iw_vt g, g2_iw_vt g2);

template <Channel_t CH> g2_iw_t chi_from_gg2(g_iw_cvt g, g2_iw_cvt g2);

g2_iw_t chi_from_gg2_PH(g_iw_vt g, g2_iw_vt g2);
g2_iw_t chi_from_gg2_PP(g_iw_vt g, g2_iw_vt g2);

} // namespace tprf
