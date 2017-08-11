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

#include <triqs/gfs.hpp>

using namespace triqs::gfs;

namespace tprf {

/// Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)
enum class Channel_t { PP, PH, PH_bar };

/// Two-particle response function type with three Matsubara frequencies
typedef gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_iw_t;

/// Short hand constant view type for g2_iw_t
typedef g2_iw_t::const_view_type g2_iw_cvt;

/// Short hand view type for g2_iw_t
typedef g2_iw_t::view_type g2_iw_vt;
  
using scalar_t = g2_iw_t::scalar_t;

/// Container type of one-particle Green and Vertex functions in Matsubara
using g_iw_t = gf<imfreq, matrix_valued>;

/// Container type of two-particle Green and Vertex functions in imaginary time
using g_tau_t = gf<imtime, matrix_valued>;

/// Container type of $\chi_3$ in Matsubara frequencies
using chi2_iw_t = gf<imfreq, tensor_valued<4>>;

/// Container type of $\chi_3$ in imaginary time
using chi2_tau_t = gf<imtime, tensor_valued<4>>;

/// Container type of $\chi_3$ in Matsubara frequencies
using chi3_iw_t = gf<cartesian_product<imfreq, imfreq>, tensor_valued<4>>;

/// Container type of $\chi_3$ in imaginary time
using chi3_tau_t = gf<cartesian_product<imtime, imtime>, tensor_valued<4>>;

/// Container type of two-particle Green and Vertex functions in Matsubara
/// frequencies
using chi4_iw_t =
    gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;

/// Container type of two-particle Green and Vertex functions in imaginary time
using chi4_tau_t =
    gf<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>;

} // namespace tprf
