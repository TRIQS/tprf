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

#include "gw.hpp"
#include "common.hpp"
#include "../mpi.hpp"

namespace triqs_tprf {

  template<typename g_out_t, typename g_in_t, typename mesh_t>  
  g_out_t dlr_on_immesh_template(g_in_t g_c, mesh_t mesh) {
    g_out_t g_out(mesh, g_c.target_shape());
    for (auto p : mesh) {
      g_out[p] = g_c(p);
    }
    return g_out;
  }

  g_w_t dlr_on_imfreq(g_Dc_cvt g_c, mesh::imfreq wmesh) {
    return dlr_on_immesh_template<g_w_t, g_Dc_cvt, mesh::imfreq>(g_c, wmesh);
  }

  chi_w_t dlr_on_imfreq(chi_Dc_cvt chi_c, mesh::imfreq wmesh) {
    return dlr_on_immesh_template<chi_w_t, chi_Dc_cvt, mesh::imfreq>(chi_c, wmesh);
  }

  g_t_t dlr_on_imtime(g_Dc_cvt g_c, mesh::imtime tmesh) {
    return dlr_on_immesh_template<g_t_t, g_Dc_cvt, mesh::imtime>(g_c, tmesh);
  }

  chi_t_t dlr_on_imtime(chi_Dc_cvt chi_c, mesh::imtime tmesh) {
    return dlr_on_immesh_template<chi_t_t, chi_Dc_cvt, mesh::imtime>(chi_c, tmesh);
  }

  std::tuple<chi_wk_t, chi_k_t> split_into_dynamic_wk_and_constant_k(chi_wk_cvt chi_wk) {

    auto _               = all_t{};
    auto wmesh = std::get<0>(chi_wk.mesh());
    auto kmesh = std::get<1>(chi_wk.mesh());

    chi_wk_t chi_dyn_wk(chi_wk.mesh(), chi_wk.target_shape());
    chi_dyn_wk() = 0.0;
    chi_k_t chi_const_k(kmesh, chi_wk.target_shape());
    chi_const_k() = 0.0;

    auto arr = mpi_view(kmesh);

#pragma omp parallel for 
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &k = arr[idx];
      
      auto chi_w = chi_wk[_, k];
      auto tail  = std::get<0>(fit_tail(chi_w));

      for (auto [a, b, c, d] : chi_wk.target_indices()) chi_const_k[k](a, b, c, d) = tail(0, a, b, c, d);

      for (auto w : wmesh) chi_dyn_wk[w, k] = chi_wk[w, k] - chi_const_k[k];
    }

    chi_dyn_wk = mpi::all_reduce(chi_dyn_wk);
    chi_const_k = mpi::all_reduce(chi_const_k);
    return {chi_dyn_wk, chi_const_k};
  }

  template<typename g_out_t, typename g_dyn_t, typename g_stat_t>  
  g_out_t add_dynamic_and_static_template(g_dyn_t g_dyn_wk, g_stat_t g_stat_k) {
  
  g_out_t g_wk(g_dyn_wk.mesh(), g_dyn_wk.target_shape());
  g_wk() = 0.0;

  auto arr = mpi_view(g_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
      auto &[w, k] = arr[idx];

      //for (const auto &[a, b] : g_wk.target_indices()) { g_wk[w, k](a, b) = g_dyn_wk[w, k](a, b) + g_stat_k[k](a, b); }
      g_wk[w, k] = g_dyn_wk[w, k] + g_stat_k[k];
  }
  g_wk = mpi::all_reduce(g_wk);
  return g_wk;
  }

  g_fk_t add_dynamic_and_static(g_fk_t g_dyn_fk, e_k_t g_stat_k) {
    return add_dynamic_and_static_template<g_fk_t, g_fk_t, e_k_t>(g_dyn_fk, g_stat_k);
  }

  chi_fk_t add_dynamic_and_static(chi_fk_t chi_dyn_fk, chi_k_t chi_stat_k) {
    return add_dynamic_and_static_template<chi_fk_t, chi_fk_t, chi_k_t>(chi_dyn_fk, chi_stat_k);
  }

  g_wk_t add_dynamic_and_static(g_wk_t g_dyn_wk, e_k_t g_stat_k) {
    return add_dynamic_and_static_template<g_wk_t, g_wk_t, e_k_t>(g_dyn_wk, g_stat_k);
  }

  chi_wk_t add_dynamic_and_static(chi_wk_t chi_dyn_wk, chi_k_t chi_stat_k) {
    return add_dynamic_and_static_template<chi_wk_t, chi_wk_t, chi_k_t>(chi_dyn_wk, chi_stat_k);
  }

  g_Dwk_t add_dynamic_and_static(g_Dwk_t g_dyn_wk, e_k_t g_stat_k) {
    return add_dynamic_and_static_template<g_Dwk_t, g_Dwk_t, e_k_t>(g_dyn_wk, g_stat_k);
  }

  chi_Dwk_t add_dynamic_and_static(chi_Dwk_t chi_dyn_wk, chi_k_t chi_stat_k) {
    return add_dynamic_and_static_template<chi_Dwk_t, chi_Dwk_t, chi_k_t>(chi_dyn_wk, chi_stat_k);
  }

//  g_fk_t add_dynamic_fk_and_static_k(g_fk_t g_dyn_fk, e_k_t g_stat_k) {
//  
//  g_fk_t g_fk(g_dyn_fk.mesh(), g_dyn_fk.target_shape());
//  g_fk() = 0.0;
//
//  auto arr = mpi_view(g_fk.mesh());
//#pragma omp parallel for
//  for (int idx = 0; idx < arr.size(); idx++) {
//      auto &[f, k] = arr[idx];
//
//      //for (const auto &[a, b] : g_fk.target_indices()) { g_fk[f, k](a, b) = g_dyn_fk[f, k](a, b) + g_stat_k[k](a, b); }
//      g_fk[f, k] = g_dyn_fk[f, k] + g_stat_k[k];
//  }
//  g_fk = mpi::all_reduce(g_fk);
//  return g_fk;
//  }

  double fermi(double e) {
    if( e < 0 ) {
      return 1. / (exp(e) + 1.);
    } else {
      double exp_me = exp(-e);
      return exp_me / (1 + exp_me);
    }
  }

  double bose(double e) {
    return 1. / (exp(e) - 1.);
  }
} // namespace triqs_tprf
