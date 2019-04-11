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
  W_wk *= 0.;
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

    auto V = to_matrix(V_k[k]);
    auto PI = to_matrix(PI_wk[w, k]);

    // auto W = to_matrix(W_wk[w, k]);

    array<scalar_t, 4> W_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};
    auto W = make_matrix_view(group_indices_view(W_arr, {0, 1}, {3, 2}));

    W = V * inverse(I - PI * V);

    W_wk[w, k] = W_arr;
  }

  W_wk = mpi_all_reduce(W_wk);
  return W_wk;
}

g_wk_t gw_self_energy(chi_wk_vt W_wk, g_wk_vt g_wk) {

  auto g_tr = make_gf_from_fourier<0, 1>(g_wk);
  auto W_tr = make_gf_from_fourier<0, 1>(W_wk);

  auto sigma_tr = make_gf(g_tr);
  sigma_tr *= 0.;

  for (const auto &[t, r] : g_tr.mesh()) {
    for (const auto &[a, b, c, d] : W_tr.target_indices()) {
      sigma_tr[t, r](a, b) += W_tr[t, r](a, b, c, d) * g_tr[t, r](c, d);
    }
  }

  auto sigma_wk = make_gf_from_fourier<0, 1>(sigma_tr);

  return sigma_wk;
}

} // namespace tprf
