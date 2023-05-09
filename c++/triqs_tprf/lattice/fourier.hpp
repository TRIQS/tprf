/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation and S. Käser
 * Authors: S. Käser
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
#pragma once

#include "../types.hpp"
#include "../fourier/fourier.hpp"
#include <omp.h>
#include "../mpi.hpp"

namespace triqs_tprf {


  namespace {
    using namespace fourier;
  }

template <typename Gf_type>
auto fourier_Dwr_to_Dtr_general_target(Gf_type g_wr) {

  auto _ = all_t{};
  // Get rid of structured binding declarations in this file due to issue #11
  //auto [wmesh, rmesh] = g_wr.mesh();
  auto wmesh = std::get<0>(g_wr.mesh());
  auto rmesh = std::get<1>(g_wr.mesh());

  //auto tmesh = make_adjoint_mesh(wmesh, n_tau);
  auto tmesh = triqs::mesh::dlr_imtime(wmesh);
  auto g_tr = make_gf<prod<dlr_imtime, cyclat>>({tmesh, rmesh}, g_wr.target());

  auto r_arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < r_arr.size(); idx++) {
    auto &r = r_arr(idx);
    auto g_c = dlr_coeffs_from_dlr_imfreq(g_wr[_, r]);
    g_tr[_, r] = dlr_imtime_from_dlr_coeffs(g_c);
  }
  g_tr = mpi::all_reduce(g_tr);
  return g_tr;
}

template <typename Gf_type>
auto fourier_Dtr_to_Dwr_general_target(Gf_type g_tr) {
  
  auto _ = all_t{};
  //auto [tmesh, rmesh] = g_tr.mesh();
  auto tmesh = std::get<0>(g_tr.mesh());
  auto rmesh = std::get<1>(g_tr.mesh());

  //auto wmesh = make_adjoint_mesh(tmesh, n_w);
  auto wmesh = triqs::mesh::dlr_imfreq(tmesh);
  auto g_wr = make_gf<prod<dlr_imfreq, cyclat>>({wmesh, rmesh}, g_tr.target());

  auto r_arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < r_arr.size(); idx++) {
    auto &r = r_arr(idx);
    auto g_c = dlr_coeffs_from_dlr_imtime(g_tr[_, r]);
    g_wr[_, r] = dlr_imfreq_from_dlr_coeffs(g_c);
  }
  g_wr = mpi::all_reduce(g_wr);
  return g_wr;
}
  
template <typename Gf_type>
auto fourier_wr_to_tr_general_target(Gf_type g_wr, int n_tau = -1) {

  auto _ = all_t{};
  // Get rid of structured binding declarations in this file due to issue #11
  //auto [wmesh, rmesh] = g_wr.mesh();
  auto wmesh = std::get<0>(g_wr.mesh());
  auto rmesh = std::get<1>(g_wr.mesh());

  auto tmesh = make_adjoint_mesh(wmesh, n_tau);
  auto g_tr = make_gf<prod<imtime, cyclat>>({tmesh, rmesh}, g_wr.target());

  auto r0 = *rmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wr[_, r0]), gf_view(g_tr[_, r0]));

  auto r_arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < r_arr.size(); idx++) {
    auto &r = r_arr(idx);

    auto g_w = make_gf<imfreq>(wmesh, g_wr.target());
    auto g_t = make_gf<imtime>(tmesh, g_tr.target());
  
    g_w = g_wr[_, r];

    _fourier_with_plan<0>(gf_const_view(g_w), gf_view(g_t), p);

    g_tr[_, r] = g_t;
  }
  g_tr = mpi::all_reduce(g_tr);
  return g_tr;
}

template <typename Gf_type>
auto fourier_tr_to_wr_general_target(Gf_type g_tr, int n_w = -1) {
  
  auto _ = all_t{};
  //auto [tmesh, rmesh] = g_tr.mesh();
  auto tmesh = std::get<0>(g_tr.mesh());
  auto rmesh = std::get<1>(g_tr.mesh());

  auto wmesh = make_adjoint_mesh(tmesh, n_w);
  auto g_wr = make_gf<prod<imfreq, cyclat>>({wmesh, rmesh}, g_tr.target());

  auto r0 = *rmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_tr[_, r0]), gf_view(g_wr[_, r0]));

  auto r_arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < r_arr.size(); idx++) {
    auto &r = r_arr(idx);

    auto g_t = make_gf<imtime>(tmesh, g_tr.target());
    auto g_w = make_gf<imfreq>(wmesh, g_wr.target());
  
    g_t = g_tr[_, r];

    _fourier_with_plan<0>(gf_const_view(g_t), gf_view(g_w), p);

    g_wr[_, r] = g_w;
  }
  g_wr = mpi::all_reduce(g_wr);
  return g_wr;
}

template <typename Gf_type>
auto fourier_wk_to_wr_general_target(Gf_type g_wk) {
  
  auto _ = all_t{};

  //auto [wmesh, kmesh] = g_wk.mesh();
  auto wmesh = std::get<0>(g_wk.mesh());
  auto kmesh = std::get<1>(g_wk.mesh());

  auto rmesh = make_adjoint_mesh(kmesh);
  auto g_wr = make_gf<prod<decltype(wmesh), cyclat>>({wmesh, rmesh}, g_wk.target());

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wk[w0, _]), gf_view(g_wr[w0, _]));

  auto w_arr = mpi_view(wmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < w_arr.size(); idx++) {
    auto &w = w_arr(idx);

    auto g_k = make_gf<brzone>(kmesh, g_wk.target());
    auto g_r = make_gf<cyclat>(rmesh, g_wr.target());
  
    g_k = g_wk[w, _];

    _fourier_with_plan<0>(gf_const_view(g_k), gf_view(g_r), p);

    g_wr[w, _] = g_r;
  }
  g_wr = mpi::all_reduce(g_wr);
  return g_wr;
}

template <typename Gf_type>
auto fourier_wr_to_wk_general_target(Gf_type g_wr) {
  
  auto _ = all_t{};

  //auto [wmesh, rmesh] = g_wr.mesh();
  auto wmesh = std::get<0>(g_wr.mesh());
  auto rmesh = std::get<1>(g_wr.mesh());

  auto kmesh = make_adjoint_mesh(rmesh);
  auto g_wk = make_gf<prod<decltype(wmesh), brzone>>({wmesh, kmesh}, g_wr.target());

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wr[w0, _]), gf_view(g_wk[w0, _]));

  auto w_arr = mpi_view(wmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < w_arr.size(); idx++) {
    auto &w = w_arr(idx);

    auto g_r = make_gf<cyclat>(rmesh, g_wr.target());
    auto g_k = make_gf<brzone>(kmesh, g_wk.target());
  
    g_r = g_wr[w, _];

    _fourier_with_plan<0>(gf_const_view(g_r), gf_view(g_k), p);

    g_wk[w, _] = g_k;
  }
  g_wk = mpi::all_reduce(g_wk);
  return g_wk;
}

} // namespace triqs_tprf
