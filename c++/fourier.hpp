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

/*

// can not do: using chi4_tau_view_t = gf_view<cartesian_product<imtime, imtime,
// imtime>, tensor_valued<4>>
// see issue: https://github.com/TRIQS/triqs/issues/421
// fixme!

// Here we resort to using defines... not nice
#define chi4_tau_view_t                                                        \
  gf_view<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>
#define chi3_tau_view_t                                                        \
  gf_view<cartesian_product<imtime, imtime>, tensor_valued<4>>
#define chi2_tau_view_t gf_view<imtime, tensor_valued<4>>
#define g_tau_view_t gf_view<imtime, matrix_valued>

#define chi4_iw_view_t                                                         \
  gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>
#define chi3_iw_view_t                                                         \
  gf_view<cartesian_product<imfreq, imfreq>, tensor_valued<4>>
#define chi2_iw_view_t gf_view<imfreq, tensor_valued<4>>
#define g_iw_view_t gf_view<imfreq, matrix_valued>
*/

namespace tprf {
namespace fourier {

  /*

chi4_iw_t chi4_iw_from_tau(chi4_tau_view_t chi4_tau, int n_iw);
chi3_iw_t chi3_iw_from_tau(chi3_tau_view_t chi3_tau, int n_iw);
chi2_iw_t chi2_iw_from_tau(chi2_tau_view_t chi2_tau, int n_iw);
g_iw_t g_iw_from_tau(g_tau_view_t g_tau, int n_iw);

chi4_tau_t chi4_tau_from_iw(chi4_iw_view_t chi4_iw, int n_tau);
chi3_tau_t chi3_tau_from_iw(chi3_iw_view_t chi3_iw, int n_tau);
chi2_tau_t chi2_tau_from_iw(chi2_iw_view_t chi2_iw, int n_tau);
g_tau_t g_tau_from_iw(g_iw_view_t g_iw, int n_tau);

  */  

} // namespace fourier
} // namespace tprf

