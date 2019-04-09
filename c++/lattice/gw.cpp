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

  auto to_matrix = [](auto t) {
    array<scalar_t, 4> t_arr{t, memory_layout_t<4>{0, 1, 2, 3}};
    auto t_mat = make_matrix_view(group_indices_view(t_arr, {0, 1}, {3, 2}));
    return t_mat;
    };

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);  

    auto PI = to_matrix(PI_wk[w, k]);
    auto W = to_matrix(W_wk[w, k]);
    auto V = to_matrix(V_k[k]);

    W = V * inverse(I - PI * V);
  }

  W_wk = mpi_all_reduce(W_wk);
  return W_wk;
}

} // namespace tprf
