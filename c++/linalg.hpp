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

#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>

using namespace triqs::gfs;
using namespace triqs::arrays;

#include "types.hpp"
#include "channel_grouping.hpp"

namespace tprf {

 /// [G]^{-1}, Two-particle response-function inversion in the PH channel
  template <Channel_t CH> g2_iw_t inverse(g2_iw_cvt g);

 /// C = A * B, two-particle response-function product the PH channel
 template <Channel_t CH> g2_iw_t product(g2_iw_cvt A, g2_iw_cvt B);

 /// 1, identity two-particle response-function the PH channel
 template <Channel_t CH> g2_iw_t identity(g2_iw_cvt g);

} // namespace tprf
