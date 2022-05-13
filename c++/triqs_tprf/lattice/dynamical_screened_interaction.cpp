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

#include <nda/nda.hpp>
#include <nda/linalg/eigenelements.hpp>

#include "dynamical_screened_interaction.hpp"
#include "common.hpp"
#include "../mpi.hpp"

namespace triqs_tprf {

chi_wk_t dynamical_screened_interaction_W(chi_wk_cvt PI_wk, chi_k_cvt V_k) {

  if( std::get<1>(PI_wk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_wk = make_gf(PI_wk);
  W_wk *= 0.;
  size_t nb = PI_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_k[k]};
    array<scalar_t, 4> PI_arr{PI_wk[w, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V  = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto PI = make_matrix_view(group_indices_view(PI_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W  = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * inverse(I - PI * V);

    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

chi_fk_t dynamical_screened_interaction_W(chi_fk_cvt PI_fk, chi_k_cvt V_k) {

  if( std::get<1>(PI_fk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_fk = make_gf(PI_fk);
  W_fk *= 0.;
  size_t nb = PI_fk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_fk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[f, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_k[k]};
    array<scalar_t, 4> PI_arr{PI_fk[f, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V  = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto PI = make_matrix_view(group_indices_view(PI_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W  = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * inverse(I - PI * V);

    W_fk[f, k] = W_arr;
  }

  W_fk = mpi::all_reduce(W_fk);
  return W_fk;
}

chi_wk_t dynamical_screened_interaction_W(chi_wk_cvt PI_wk, chi_wk_cvt V_wk) {

  if( std::get<0>(PI_wk.mesh()) != std::get<0>(V_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: Matsubara meshes are not the same\n";

  if( std::get<1>(PI_wk.mesh()) != std::get<1>(V_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_wk = make_gf(PI_wk);
  W_wk *= 0.;
  size_t nb = PI_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_wk[w, k]};
    array<scalar_t, 4> PI_arr{PI_wk[w, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V  = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto PI = make_matrix_view(group_indices_view(PI_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W  = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * inverse(I - PI * V);

    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

chi_fk_t dynamical_screened_interaction_W(chi_fk_cvt PI_fk, chi_fk_cvt V_fk) {

  if( std::get<0>(PI_fk.mesh()) != std::get<0>(V_fk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: real-frequency meshes are not the same\n";

  if( std::get<1>(PI_fk.mesh()) != std::get<1>(V_fk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_fk = make_gf(PI_fk);
  W_fk *= 0.;
  size_t nb = PI_fk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_fk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[f, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_fk[f, k]};
    array<scalar_t, 4> PI_arr{PI_fk[f, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V  = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto PI = make_matrix_view(group_indices_view(PI_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W  = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * inverse(I - PI * V);

    W_fk[f, k] = W_arr;
  }

  W_fk = mpi::all_reduce(W_fk);
  return W_fk;
}





chi_wk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k) {

  if( std::get<1>(chi_wk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_wk = make_gf(chi_wk);
  W_wk *= 0.;
  size_t nb = W_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_k[k]};
    array<scalar_t, 4> chi_arr{chi_wk[w, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V   = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto chi = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W   = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * chi * V + V;
    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

chi_fk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk_cvt chi_fk, chi_k_cvt V_k) {

  if( std::get<1>(chi_fk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_fk = make_gf(chi_fk);
  W_fk *= 0.;
  size_t nb = W_fk.target_shape()[0];

  using scalar_t = chi_fk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_fk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_k[k]};
    array<scalar_t, 4> chi_arr{chi_fk[w, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V   = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto chi = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W   = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * chi * V + V;
    W_fk[w, k] = W_arr;
  }

  W_fk = mpi::all_reduce(W_fk);
  return W_fk;
}

chi_wk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_wk_cvt V_wk) {

  if( std::get<0>(chi_wk.mesh()) != std::get<0>(V_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: Matsubara meshes are not the same\n";

  if( std::get<1>(chi_wk.mesh()) != std::get<1>(V_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_wk = make_gf(chi_wk);
  W_wk *= 0.;
  size_t nb = W_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_wk[w, k]};
    array<scalar_t, 4> chi_arr{chi_wk[w, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V   = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto chi = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W   = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * chi * V + V;
    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

chi_fk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk_cvt chi_fk, chi_fk_cvt V_fk) {

  if( std::get<0>(chi_fk.mesh()) != std::get<0>(V_fk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: Real-frequency meshes are not the same\n";

  if( std::get<1>(chi_fk.mesh()) != std::get<1>(V_fk.mesh()) )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";
  
  auto W_fk = make_gf(chi_fk);
  W_fk *= 0.;
  size_t nb = W_fk.target_shape()[0];

  using scalar_t = chi_fk_t::scalar_t;
  auto I = nda::eye<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_fk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_fk[w, k]};
    array<scalar_t, 4> chi_arr{chi_fk[w, k]};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb};

    auto V   = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto chi = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));
    auto W   = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

    W = V * chi * V + V;
    W_fk[w, k] = W_arr;
  }

  W_fk = mpi::all_reduce(W_fk);
  return W_fk;
}


} // namespace triqs_tprf
