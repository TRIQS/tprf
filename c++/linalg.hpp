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

#include "channel_grouping.hpp"
#include "types.hpp"

namespace tprf {

/** Two-particle response-function inversion $[g]^{-1}$.
 
 The two-particle response function $g_{abcd}(\omega, \nu, \nu')$ 
 is cast to matrix form and inverted

 .. math::
   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

 where the mapping of target-space indices $\{a, b, c, d \}$ to $\{\alpha, \beta\}, \{\gamma, \delta\}$ is channel dependent.
 
 Storage is allocated and the inverse is returned by value.
 
 @tparam CH selects the two-particle channel
 @param g two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`
 @return :math:`[g]^{-1}` in the given channel
 @include tprf/linalg.hpp
 @note Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation
 
 */
template <Channel_t CH> g2_iw_t inverse(g2_iw_cvt g);

/** Two-particle response-function product $A * B$
 
 The two-particle response functions $A \equiv A_{abcd}(\omega, \nu, \nu')$ 
 and $B \equiv A_{abcd}(\omega, \nu, \nu')$ are cast to matrix form and their
 product is computed

 .. math::
   (A * B)_{abcd}(\omega, \nu, \nu') 
   = \sum_{\nu_1 A_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

 where the mapping of target-space indices $\{a, b, c, d \}$ to $\{\alpha, \beta\}, \{\gamma, \delta\}$ is channel dependent.

 FIXME
 
 Storage is allocated and the inverse is returned by value.
 
 @tparam CH selects the two-particle channel
 @param g two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`
 @return :math:`[g]^{-1}` in the given channel
 @include tprf/linalg.hpp
 @note Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation
 
 */
template <Channel_t CH> g2_iw_t product(g2_iw_cvt A, g2_iw_cvt B);

/// 1, identity two-particle response-function the PH channel
template <Channel_t CH> g2_iw_t identity(g2_iw_cvt g);

} // namespace tprf
