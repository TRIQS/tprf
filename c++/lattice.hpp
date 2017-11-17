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

#include <triqs/lattice/brillouin_zone.hpp>
#include <triqs/lattice/tight_binding.hpp>
#include <triqs/utility/mini_vector.hpp>

namespace tprf {

gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh);
gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma);

gr_iw_t gr_from_gk(gk_iw_vt gk);
gk_iw_t gk_from_gr(gr_iw_vt gr, brillouin_zone bz);

chi0r_t chi0r_from_gr_PH(int nw, int nnu, gr_iw_vt gr);

chi0r_t chi0r_from_chi0q(chi0q_vt chi0q);
chi0q_t chi0q_from_chi0r(chi0r_vt chi0r, brillouin_zone bz);

gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu(chi0q_t chi0q);

chiq_t chiq_from_chi0q_and_gamma_PH(chi0q_vt chi0q, g2_iw_vt gamma_ph);
  
} // namespace tprf
