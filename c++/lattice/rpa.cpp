/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018, The Simons Foundation
 * Author: Hugo U. R. Strand
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

#include "common.hpp"
#include "rpa.hpp"

namespace tprf {

  chi_wk_t solve_rpa_PH(chi_wk_vt chi0_wk, array_view<std::complex<double>, 4> U_arr) {

    using scalar_t = chi_wk_t::scalar_t;

    size_t nb = chi0_wk.target_shape()[0];
    auto wmesh = std::get<0>(chi0_wk.mesh());
    auto kmesh = std::get<1>(chi0_wk.mesh());

    chi_wk_t chi_wk{{wmesh, kmesh}, chi0_wk.target_shape()};

    auto U = make_matrix_view(group_indices_view(U_arr, {0, 1}, {3, 2})); // PH grouping, from c+c+cc
    auto I = make_unit_matrix<scalar_t>(U.shape()[0]);
    
    for (auto const &w : wmesh) {
      for (auto const &k : kmesh) {

	array<scalar_t, 4> chi_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};
	array<scalar_t, 4> chi0_arr{chi0_wk[w, k], memory_layout_t<4>{0, 1, 2, 3}};
	
	auto chi = make_matrix_view(group_indices_view(chi_arr, {0, 1}, {3, 2})); // PH grouping FIXME!
	auto chi0 = make_matrix_view(group_indices_view(chi0_arr, {0, 1}, {3, 2})); // PH grouping FIXME!

	chi = inverse(I - chi0 * U) * chi0; // Inverted BSE specialized for rpa

	chi_wk[w, k] = chi_arr; // assign back using the array_view
	
      }
    }

    return chi_wk;
    
  }
  
}
