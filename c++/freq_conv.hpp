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

g_iw_t block_iw_AB_to_matrix_valued(b_g_iw_vt bg_AB);
  
void block_3nu_AABB_to_tensor_valued(b_g2_iw_vt bg2_AABB, g2_iw_vt g2);

void get_magnetic_component(g2_iw_vt g2, g2_iw_vt g2_m);

template <Channel_t CH> void from_3nu(g2_iw_vt g2_ch, g2_iw_cvt g2);

void from_3nu_PH(g2_iw_vt g2_ch, g2_iw_vt g2);
void from_3nu_PH_bar(g2_iw_vt g2_ch, g2_iw_vt g2);
void from_3nu_PP(g2_iw_vt g2_ch, g2_iw_vt g2);

} // namespace tprf
