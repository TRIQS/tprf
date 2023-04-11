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
#include <omp.h>
#include "../mpi.hpp"

namespace triqs_tprf {

  template<typename CHI_T, typename CHI_VT>
  CHI_T solve_rpa_PH(CHI_VT chi0_wk, array_contiguous_view<std::complex<double>, 4> U_arr) {

    using scalar_t = chi_wk_t::scalar_t;

    size_t nb = chi0_wk.target_shape()[0];

    auto chi_wk = make_gf(chi0_wk);
    chi_wk *= 0;

    // PH grouping of the vertex, from cc+cc+, permuting the last two indices.
    auto U = make_matrix_view(group_indices_view(U_arr, idx_group<0, 1>, idx_group<3, 2>));

    auto I = nda::eye<scalar_t>(U.shape()[0]);

    auto meshes_mpi = mpi_view(chi0_wk.mesh());

#pragma omp parallel for
  for (unsigned int idx = 0; idx < meshes_mpi.size(); idx++){
    auto &[w, k] = meshes_mpi[idx];

    array<scalar_t, 4> chi_arr{nb, nb, nb, nb};
    array<scalar_t, 4> chi0_arr{chi0_wk[w, k]};

    // PH grouping (permuting last two indices)
    auto chi  = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto chi0 = make_matrix_view(group_indices_view(chi0_arr, idx_group<0, 1>, idx_group<3, 2>));

    chi = inverse(I - chi0 * U) * chi0; // Inverted BSE specialized for rpa

    chi_wk[w, k] = chi_arr;             // assign back using the array_view
    }
  chi_wk = mpi::all_reduce(chi_wk);

  return chi_wk;
  }

  chi_wk_t solve_rpa_PH(chi_wk_vt chi0_wk, array_contiguous_view<std::complex<double>, 4> U_arr) {
    return solve_rpa_PH<chi_wk_t, chi_wk_vt>(chi0_wk, U_arr);
  }
      
  chi_fk_t solve_rpa_PH(chi_fk_vt chi0_fk, array_contiguous_view<std::complex<double>, 4> U_arr) {
    return solve_rpa_PH<chi_fk_t, chi_fk_vt>(chi0_fk, U_arr);
  }

} // namespace triqs_tprf
