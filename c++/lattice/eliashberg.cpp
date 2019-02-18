/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation and S. Käser
 * Authors: H. U.R. Strand, S. Käser
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

#include "eliashberg.hpp"

namespace tprf {

// Helper function computing F = GG \Delta
  
gk_iw_t eliashberg_g_delta_g_product(gk_iw_vt g_wk, gk_iw_vt delta_wk) {

  auto [wmesh, kmesh] = delta_wk.mesh();

  auto F_wk = make_gf(delta_wk);
  F_wk *= 0.;

  for (const auto [w, k] : delta_wk.mesh())
    for (auto [A, B] : F_wk.target_indices())
      for (auto [c, d] : delta_wk.target_indices())
        F_wk[w, k](A, B) +=
            g_wk[w, k](A, c) * g_wk[-w, -k](B, d) * delta_wk[w, k](c, d);

  return F_wk;
}

gk_iw_t eliashberg_product(chi_wk_vt Gamma_pp, gk_iw_vt g_wk,
                           gk_iw_vt delta_wk) {

  auto [wmesh, kmesh] = delta_wk.mesh();

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);
  
  auto delta_wk_out = make_gf(delta_wk);
  delta_wk_out *= 0.;

  for (const auto [w, k] : delta_wk.mesh())
    for (const auto [n, q] : delta_wk.mesh())
      for (auto [A, a, b, B] : Gamma_pp.target_indices())
        delta_wk_out[w, k](a, b) +=
            Gamma_pp[w - n, k - q](A, a, b, B) * F_wk[n, k](A, B);

  delta_wk_out /= -wmesh.domain().beta;

  return delta_wk_out;
}

} // namespace tprf
