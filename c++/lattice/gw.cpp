/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation
 * Authors: H. U.R. Strand
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

#include "gw.hpp"
#include "../mpi.hpp"

namespace tprf {

chi_wk_t screened_interaction_W(chi_wk_vt PI_wk, chi_k_vt V_k) {

  auto W_wk = make_gf(PI_wk);
  size_t nb = PI_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = make_unit_matrix<scalar_t>(nb * nb);

  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);  

  //for (auto const &[w, k] : W_wk.mesh()) {
    
    array<scalar_t, 4> PI_arr{PI_wk[w, k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> V_arr{V_k[k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};

    // PH grouping (permuting last two indices)
    auto PI = make_matrix_view(group_indices_view(PI_arr, {0, 1}, {3, 2}));
    auto V = make_matrix_view(group_indices_view(V_arr, {0, 1}, {3, 2}));
    auto W = make_matrix_view(group_indices_view(W_arr, {0, 1}, {3, 2}));

    W = V * inverse(I - PI * V);

    W_wk[w, k] = W_arr; // assign back using the array_view
  }
  
  return W_wk;
}

} // namespace tprf
