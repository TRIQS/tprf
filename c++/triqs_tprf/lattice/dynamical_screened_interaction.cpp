/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, The Simons Foundation
 * Authors: H. U.R. Strand, Y. in 't Veld, N. Wentzell
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

  enum SusceptibilityType { bubble, generalized };
 
  template <SusceptibilityType susType, typename chi_t, typename v_t> auto screened_interaction_from_generic_susceptibility(chi_t &chi, v_t &V) {

    auto const &[freqmesh, kmesh] = chi.mesh();

    auto const &vkmesh = [&V]() -> auto & {
      if constexpr (v_t::arity == 1)
        return V.mesh();
      else
        return std::get<1>(V.mesh());
    }();

    if (kmesh != vkmesh) TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";

    auto W    = make_gf(chi);
    W()       = 0.;
    size_t nb = chi.target_shape()[0];

    using scalar_t = typename chi_t::scalar_t;
    auto I         = nda::eye<scalar_t>(nb * nb);

    // MPI and openMP parallell loop
    auto arr = mpi_view(W.mesh());
#pragma omp parallel for
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &[w, k] = arr[idx];

      array<scalar_t, 4> V_arr;
      if constexpr (v_t::arity == 1)
        V_arr = V[k];
      else
        V_arr = V[w, k];

      array<scalar_t, 4> chi_arr{chi[w, k]};
      array<scalar_t, 4> W_arr{nb, nb, nb, nb};

      auto V_mat   = make_matrix_view(group_indices_view(V_arr, idx_group<0, 1>, idx_group<3, 2>));
      auto chi_mat = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));
      auto W_mat   = make_matrix_view(group_indices_view(W_arr, idx_group<0, 1>, idx_group<3, 2>));

      if constexpr (susType == false)
        W_mat = V_mat * inverse(I - chi_mat * V_mat);
      else
        W_mat = V_mat * chi_mat * V_mat + V_mat;

      W[w, k] = W_arr;
    }

    W = mpi::all_reduce(W);
    return W;
  }

  template <typename chi_t>
    auto dynamical_screened_interaction_W_opt_bubble(chi_t chi_wk, chi_k_cvt V_k) {
    
    auto const &[wmesh, kmesh] = chi_wk.mesh();
    auto const &v_kmesh = V_k.mesh();

    if (kmesh != v_kmesh) TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W: k-space meshes are not the same\n";

    auto W_wk    = make_gf(chi_wk);
    W_wk()       = 0.;
    size_t nb = chi_wk.target_shape()[0];

    using scalar_t = typename chi_wk_cvt::scalar_t;
    auto I         = nda::eye<scalar_t>(nb * nb);
    
    // MPI and openMP parallell loop
    auto arr = mpi_view(W_wk.mesh());
    
#pragma omp parallel for 
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &[w, k] = arr[idx];

      array<scalar_t, 4> denom{nb, nb, nb, nb};
      array<scalar_t, 4> inv_denom{nb, nb, nb, nb};
      
      auto denom_mat = make_matrix_view(group_indices_view(denom, idx_group<0, 1>, idx_group<3, 2>));
      auto inv_denom_mat = make_matrix_view(group_indices_view(inv_denom, idx_group<0, 1>, idx_group<3, 2>));

      denom_mat = I;

      for (auto const &[a, b, c, d] : V_k.target_indices())
	for( auto e : range(nb) ) 
	  for( auto f : range(nb) ) 
	    denom(a,b,c,d) -= chi_wk[w,k](a,b,e,f)*V_k[k](e,f,c,d);
      
      inv_denom_mat = inverse(denom_mat);

      for (auto const &[a, b, c, d] : V_k.target_indices())
	for( auto e : range(nb) ) 
	  for( auto f : range(nb) ) 
	    W_wk[w, k](a,b,c,d) += V_k[k](a,b,e,f) * inv_denom(e,f,c,d);      
    }

    W_wk = mpi::all_reduce(W_wk);
    return W_wk;
  }

  chi_wk_t dynamical_screened_interaction_W(chi_wk_cvt PI_wk, chi_k_cvt V_k) {
    return dynamical_screened_interaction_W_opt_bubble<chi_wk_cvt>(PI_wk, V_k);
  }

  chi_Dwk_t dynamical_screened_interaction_W(chi_Dwk_cvt PI_wk, chi_k_cvt V_k) {
    return dynamical_screened_interaction_W_opt_bubble<chi_Dwk_cvt>(PI_wk, V_k);
  }
  
  chi_fk_t dynamical_screened_interaction_W(chi_fk_cvt PI_fk, chi_k_cvt V_k) {
    return screened_interaction_from_generic_susceptibility<bubble>(PI_fk, V_k);
  }

  chi_wk_t dynamical_screened_interaction_W(chi_wk_cvt PI_wk, chi_wk_cvt V_wk) {
    return screened_interaction_from_generic_susceptibility<bubble>(PI_wk, V_wk);
  }

  chi_Dwk_t dynamical_screened_interaction_W(chi_Dwk_cvt PI_wk, chi_Dwk_cvt V_wk) {
    return screened_interaction_from_generic_susceptibility<bubble>(PI_wk, V_wk);
  }

  chi_fk_t dynamical_screened_interaction_W(chi_fk_cvt PI_fk, chi_fk_cvt V_fk) {
    return screened_interaction_from_generic_susceptibility<bubble>(PI_fk, V_fk);
  }

  chi_wk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k) {
    return screened_interaction_from_generic_susceptibility<generalized>(chi_wk, V_k);
  }

  chi_Dwk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_Dwk_cvt chi_wk, chi_k_cvt V_k) {
    return screened_interaction_from_generic_susceptibility<generalized>(chi_wk, V_k);
  }

  chi_fk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk_cvt chi_fk, chi_k_cvt V_k) {
    return screened_interaction_from_generic_susceptibility<generalized>(chi_fk, V_k);
  }

  chi_wk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_wk_cvt V_wk) {
    return screened_interaction_from_generic_susceptibility<generalized>(chi_wk, V_wk);
  }

  chi_Dwk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_Dwk_cvt chi_wk, chi_Dwk_cvt V_wk) {
    return screened_interaction_from_generic_susceptibility<generalized>(chi_wk, V_wk);
  }
  
  chi_fk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk_cvt chi_fk, chi_fk_cvt V_fk) {
    return screened_interaction_from_generic_susceptibility<generalized>(chi_fk, V_fk);
  }

} // namespace triqs_tprf
