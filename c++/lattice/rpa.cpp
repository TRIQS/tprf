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

#include "rpa.hpp"
#include "common.hpp"

namespace tprf {

chi_wk_t solve_rpa_PH(chi_wk_vt chi0_wk,
                      array_view<std::complex<double>, 4> U_arr) {

  using scalar_t = chi_wk_t::scalar_t;

  size_t nb = chi0_wk.target_shape()[0];
  auto wmesh = std::get<0>(chi0_wk.mesh());
  auto kmesh = std::get<1>(chi0_wk.mesh());

  chi_wk_t chi_wk{{wmesh, kmesh}, chi0_wk.target_shape()};

  // PH grouping of the vertex, from cc+cc+, permuting the last two indices.
  auto U = make_matrix_view(group_indices_view(U_arr, {0, 1}, {3, 2}));

  auto I = make_unit_matrix<scalar_t>(U.shape()[0]);

  for (auto const &w : wmesh) {
    for (auto const &k : kmesh) {

      array<scalar_t, 4> chi_arr{nb, nb, nb, nb,
                                 memory_layout_t<4>{0, 1, 2, 3}};
      array<scalar_t, 4> chi0_arr{chi0_wk[w, k],
                                  memory_layout_t<4>{0, 1, 2, 3}};

      // PH grouping (permuting last two indices)
      auto chi = make_matrix_view(group_indices_view(chi_arr, {0, 1}, {3, 2}));
      auto chi0 =
          make_matrix_view(group_indices_view(chi0_arr, {0, 1}, {3, 2}));

      chi = inverse(I - chi0 * U) * chi0; // Inverted BSE specialized for rpa

      chi_wk[w, k] = chi_arr; // assign back using the array_view 
    }
  }

  return chi_wk;
}

chi_wk_t solve_rpa_spin(chi_wk_vt chi0_wk,
                      array_view<std::complex<double>, 4> U_arr) {

  using scalar_t = chi_wk_t::scalar_t;

  size_t nb = chi0_wk.target_shape()[0];
  auto wmesh = std::get<0>(chi0_wk.mesh());
  auto kmesh = std::get<1>(chi0_wk.mesh());

  chi_wk_t chi_wk{{wmesh, kmesh}, chi0_wk.target_shape()};

  // Normal matrix Product
  auto U = make_matrix_view(group_indices_view(U_arr, {0, 1}, {2, 3}));

  auto I = make_unit_matrix<scalar_t>(U.shape()[0]);

  for (auto const &w : wmesh) {
    for (auto const &k : kmesh) {

      array<scalar_t, 4> chi_arr{nb, nb, nb, nb,
                                 memory_layout_t<4>{0, 1, 2, 3}};
      array<scalar_t, 4> chi0_arr{chi0_wk[w, k],
                                  memory_layout_t<4>{0, 1, 2, 3}};

      auto chi = make_matrix_view(group_indices_view(chi_arr, {0, 1}, {2, 3}));
      auto chi0 =
          make_matrix_view(group_indices_view(chi0_arr, {0, 1}, {2, 3}));

      chi = inverse(I - chi0 * U) * chi0; // Inverted BSE specialized for rpa

      chi_wk[w, k] = chi_arr; // assign back using the array_view
    }
  }

  return chi_wk;
}

chi_wk_t solve_rpa_charge(chi_wk_vt chi0_wk,
                      array_view<std::complex<double>, 4> U_arr) {

  using scalar_t = chi_wk_t::scalar_t;

  size_t nb = chi0_wk.target_shape()[0];
  auto wmesh = std::get<0>(chi0_wk.mesh());
  auto kmesh = std::get<1>(chi0_wk.mesh());

  chi_wk_t chi_wk{{wmesh, kmesh}, chi0_wk.target_shape()};

  // Normal matrix product
  auto U = make_matrix_view(group_indices_view(U_arr, {0, 1}, {2, 3}));

  auto I = make_unit_matrix<scalar_t>(U.shape()[0]);

  for (auto const &w : wmesh) {
    for (auto const &k : kmesh) {

      array<scalar_t, 4> chi_arr{nb, nb, nb, nb,
                                 memory_layout_t<4>{0, 1, 2, 3}};
      array<scalar_t, 4> chi0_arr{chi0_wk[w, k],
                                  memory_layout_t<4>{0, 1, 2, 3}};

      auto chi = make_matrix_view(group_indices_view(chi_arr, {0, 1}, {2, 3}));
      auto chi0 =
          make_matrix_view(group_indices_view(chi0_arr, {0, 1}, {2, 3}));

      chi = inverse(I + chi0 * U) * chi0; // Inverted BSE specialized for rpa

      chi_wk[w, k] = chi_arr; // assign back using the array_view
    }
  }

  return chi_wk;
}

} // namespace tprf
