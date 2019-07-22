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

#include "eliashberg.hpp"
#include "common.hpp"
#include "../mpi.hpp"
#include <triqs/utility/timer.hpp>
#include "../fourier/fourier.hpp"

namespace triqs_tprf {

  namespace {
    using namespace fourier;
  }

// Helper function computing F = GG \Delta

g_wk_t eliashberg_g_delta_g_product(g_wk_vt g_wk, g_wk_vt delta_wk) {
  triqs::utility::timer t_all, t_parallel;
  t_all.start();

  size_t norb = g_wk.target_shape()[0];

  auto [wmesh, kmesh] = delta_wk.mesh();
  auto wmesh_gf = std::get<0>(g_wk.mesh());

  if (wmesh.size() > wmesh_gf.size())
      TRIQS_RUNTIME_ERROR << "The size of the Matsubara frequency mesh of the Green's function"
          " (" << wmesh_gf.size() << ") must be atleast the size of the mesh of Delta (" <<
          wmesh.size() << ").";

  auto F_wk = make_gf(delta_wk);
  F_wk *= 0.;

/* The rest of this function contains a lot boiler plate code due to issue
   #725 in the TRIQS library and not yet avaible functionality to use 
   'pragma omp parallel loop for' over mesh objects.
   It will be changed later
*/
  array<typename gf_mesh<brillouin_zone>::mesh_point_t, 1> k_arr(kmesh.size());
    auto k_iter = kmesh.begin();
    for (auto idx : range(0, kmesh.size())) {
        auto k  = *k_iter;
        k_arr(idx) = k;
        k_iter++;
    }

  auto _ = all_t{};
  t_parallel.start();
#pragma omp parallel for
 for(int idx_k = 0; idx_k < kmesh.size(); idx_k++){
   auto k = k_arr(idx_k);

   auto g_left_w = make_gf<imfreq>(wmesh, g_wk.target());
   auto g_right_w = make_gf<imfreq>(wmesh, g_wk.target());
   auto delta_w = make_gf<imfreq>(wmesh, delta_wk.target());
   auto F_w = make_gf<imfreq>(wmesh, F_wk.target());

   g_left_w = g_wk[_, k]; 
   g_right_w = g_wk[_, -k]; 
   delta_w = delta_wk[_, k];

  for (const auto w : wmesh) {
    for (auto [A, B] : F_wk.target_indices())
      for (auto [c, d] : delta_wk.target_indices())
        F_w[w](A, B) +=
        g_left_w[w](A, c) * g_right_w[-w](B, d) * delta_w[w](c, d);
     }
     F_wk[_, k] = F_w;
   }
   t_parallel.stop();
   t_all.stop();
   std::cout << "all:\t" << double(t_all) << "\n" \
              << "parallel:\t" << double(t_parallel) << "\n" \
              << "sequential:\t" << double(t_all - t_parallel) << "\n";

  return F_wk;
}

g_wk_t eliashberg_product(chi_wk_vt Gamma_pp, g_wk_vt g_wk,
                       g_wk_vt delta_wk) {

  auto [wmesh, kmesh] = delta_wk.mesh();
  auto gamma_wmesh = std::get<0>(Gamma_pp.mesh());

  if (2*wmesh.size() > gamma_wmesh.size())
      TRIQS_RUNTIME_ERROR << "The size of the Matsubara frequency mesh of Gamma"
          " (" << gamma_wmesh.size() << ") must be atleast TWICE the size of the mesh of Delta (" <<
          wmesh.size() << ").";

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);

  auto delta_wk_out = make_gf(delta_wk);
  delta_wk_out *= 0.;
    
  for (const auto [w, k] : delta_wk.mesh())
    for (const auto [n, q] : delta_wk.mesh())
      for (auto [A, a, B, b] : Gamma_pp.target_indices())
        delta_wk_out[w, k](a, b) +=
            Gamma_pp[w-n, k - q](A, a, B, b) * F_wk[n, q](A, B);

  delta_wk_out /= -(wmesh.domain().beta * kmesh.size());
  
  return delta_wk_out;
}

std::tuple<chi_wk_vt, chi_k_vt> split_into_dynamic_wk_and_constant_k(chi_wk_vt Gamma_pp) {

  auto _ = all_t{};
  auto [wmesh, kmesh] = Gamma_pp.mesh();
    
  // Fit infinite frequency value
  auto Gamma_pp_dyn_wk = make_gf(Gamma_pp);

  chi_k_vt Gamma_pp_const_k = make_gf(kmesh, Gamma_pp.target());

  for (const auto k : kmesh) {
    auto Gamma_w = Gamma_pp[_, k];
    auto [tail, err] = fit_tail(Gamma_w);
    for (auto [a, b, c, d] : Gamma_pp.target_indices())
      Gamma_pp_const_k[k](a, b, c, d) = tail(0, a, b, c, d);
    for( const auto w : wmesh ) Gamma_pp_dyn_wk[w, k] = Gamma_pp[w, k] - Gamma_pp_const_k[k];
  }

    return {Gamma_pp_dyn_wk, Gamma_pp_const_k};
}

std::tuple<chi_tr_vt, chi_r_vt> dynamic_and_constant_to_tr(chi_wk_vt Gamma_pp_dyn_wk, 
                                                            chi_k_vt Gamma_pp_const_k) {

    auto Gamma_pp_dyn_tr = make_gf_from_fourier<0, 1>(Gamma_pp_dyn_wk);
    auto Gamma_pp_const_r = make_gf_from_fourier<0>(Gamma_pp_const_k);

    return {Gamma_pp_dyn_tr, Gamma_pp_const_r};
}

e_k_t eliashberg_constant_gamma_f_product(chi_r_vt Gamma_pp_const_r, g_tr_t F_tr) {

  auto _ = all_t{};

  auto delta_r_out = make_gf(std::get<1>(F_tr.mesh()), F_tr.target());
  delta_r_out *= 0.;

  for (const auto r : std::get<1>(F_tr.mesh())) {
    auto F_t = F_tr[_, r];
    for (auto [A, a, B, b] : Gamma_pp_const_r.target_indices())
        delta_r_out[r](a, b) += -Gamma_pp_const_r[r](A, a, B, b) * F_t(0)(A, B);
  }

  auto delta_k_out = make_gf_from_fourier<0>(delta_r_out);

  return delta_k_out;
}

// BOILER PLATE CODE FOR FFT STARTS HERE ============================================================ 

g_tr_t g_tr_from_g_wr(g_wr_cvt g_wr, int nt = -1) {
  
  auto _ = all_t{};
  int nb = g_wr.target().shape()[0];

  auto [wmesh, rmesh] = g_wr.mesh();

  auto tmesh = make_adjoint_mesh(wmesh, nt);
  g_tr_t g_tr{{tmesh, rmesh}, {nb, nb}};

  auto r0 = *rmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wr[_, r0]), gf_view(g_tr[_, r0]));

  array<typename gf_mesh<cyclic_lattice>::mesh_point_t, 1> r_arr(rmesh.size());
    auto r_iter = rmesh.begin();
    for (auto idx : range(0, rmesh.size())) {
        auto r  = *r_iter;
        r_arr(idx) = r;
        r_iter++;
    }

#pragma omp parallel for 
  for (int idx = 0; idx < r_arr.size(); idx++) {
    auto &r = r_arr(idx);

    auto g_w = make_gf<imfreq>(wmesh, g_wr.target());
    auto g_t = make_gf<imtime>(tmesh, g_tr.target());
  
    g_w = g_wr[_, r];

    _fourier_with_plan<0>(gf_const_view(g_w), gf_view(g_t), p);

    g_tr[_, r] = g_t;
  }
  return g_tr;
}

g_wr_t g_wr_from_g_wk(g_wk_cvt g_wk) {
  
  auto _ = all_t{};
  int nb = g_wk.target().shape()[0];

  auto [wmesh, kmesh] = g_wk.mesh();

  auto rmesh = make_adjoint_mesh(kmesh);
  g_wr_t g_wr{{wmesh, rmesh}, {nb, nb}};

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wk[w0, _]), gf_view(g_wr[w0, _]));

  array<typename gf_mesh<imfreq>::mesh_point_t, 1> w_arr(wmesh.size());
    auto w_iter = wmesh.begin();
    for (auto idx : range(0, wmesh.size())) {
        auto w  = *w_iter;
        w_arr(idx) = w;
        w_iter++;
    }

#pragma omp parallel for 
  for (int idx = 0; idx < w_arr.size(); idx++) {
    auto &w = w_arr(idx);

    auto g_k = make_gf<brillouin_zone>(kmesh, g_wk.target());
    auto g_r = make_gf<cyclic_lattice>(rmesh, g_wr.target());
  
    g_k = g_wk[w, _];

    _fourier_with_plan<0>(gf_const_view(g_k), gf_view(g_r), p);

    g_wr[w, _] = g_r;
  }
  return g_wr;
}

g_wr_t g_wr_from_g_tr(g_tr_cvt g_tr) {
  
  auto _ = all_t{};
  int nb = g_tr.target().shape()[0];

  auto [tmesh, rmesh] = g_tr.mesh();

  auto wmesh = make_adjoint_mesh(tmesh);
  g_wr_t g_wr{{wmesh, rmesh}, {nb, nb}};

  auto r0 = *rmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_tr[_, r0]), gf_view(g_wr[_, r0]));

  array<typename gf_mesh<cyclic_lattice>::mesh_point_t, 1> r_arr(rmesh.size());
    auto r_iter = rmesh.begin();
    for (auto idx : range(0, rmesh.size())) {
        auto r  = *r_iter;
        r_arr(idx) = r;
        r_iter++;
    }

#pragma omp parallel for 
  for (int idx = 0; idx < r_arr.size(); idx++) {
    auto &r = r_arr(idx);

    auto g_t = make_gf<imtime>(tmesh, g_tr.target());
    auto g_w = make_gf<imfreq>(wmesh, g_wr.target());
  
    g_t = g_tr[_, r];

    _fourier_with_plan<0>(gf_const_view(g_t), gf_view(g_w), p);

    g_wr[_, r] = g_w;
  }
  return g_wr;
}

g_wk_t g_wk_from_g_wr(g_wr_cvt g_wr) {
  
  auto _ = all_t{};
  int nb = g_wr.target().shape()[0];

  auto [wmesh, rmesh] = g_wr.mesh();

  auto kmesh = make_adjoint_mesh(rmesh);
  g_wk_t g_wk{{wmesh, kmesh}, {nb, nb}};

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wr[w0, _]), gf_view(g_wk[w0, _]));

  array<typename gf_mesh<imfreq>::mesh_point_t, 1> w_arr(wmesh.size());
    auto w_iter = wmesh.begin();
    for (auto idx : range(0, wmesh.size())) {
        auto w  = *w_iter;
        w_arr(idx) = w;
        w_iter++;
    }

#pragma omp parallel for 
  for (int idx = 0; idx < w_arr.size(); idx++) {
    auto &w = w_arr(idx);

    auto g_r = make_gf<cyclic_lattice>(rmesh, g_wr.target());
    auto g_k = make_gf<brillouin_zone>(kmesh, g_wk.target());
  
    g_r = g_wr[w, _];

    _fourier_with_plan<0>(gf_const_view(g_r), gf_view(g_k), p);

    g_wk[w, _] = g_k;
  }
  return g_wk;
}

// BOILER PLATE CODE FOR FFT ENDS HERE ============================================================ 

g_wk_t eliashberg_product_fft(chi_tr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r,
                                g_wk_vt g_wk, g_wk_vt delta_wk) {

  triqs::utility::timer t_g_delta_g_product, t_fft_F, t_dynamic_product, t_fft_delta, t_constant_product, t_combine;
    
  t_g_delta_g_product.start();
  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);
  t_g_delta_g_product.stop();

  t_fft_F.start();
  //auto F_tr = make_gf_from_fourier<0, 1>(F_wk);
  auto F_wr = g_wr_from_g_wk(F_wk);
  auto F_tr = g_tr_from_g_wr(F_wr);
  t_fft_F.stop();

  auto [tmesh, rmesh] = F_tr.mesh();

  // Dynamic part
  auto delta_tr_out = make_gf(F_tr);
  delta_tr_out *= 0.;

  auto tmesh_gamma = std::get<0>(Gamma_pp_dyn_tr.mesh());

  // Test if the tau meshs of delta and gamma are compatible. If not raise an error, because
  // it would lead to wrong results.
  if (tmesh.size() != tmesh_gamma.size()) 
      TRIQS_RUNTIME_ERROR << "The size of the imaginary time mesh of Gamma"
          " (" << tmesh_gamma.size() << ") must be the size of the mesh of Delta (" <<
          tmesh.size() << ").";

  t_dynamic_product.start();
/* This function contains a lot boiler plate code due to issue
   #725 in the TRIQS library and not yet avaible functionality to use 
   'pragma omp parallel loop for' over mesh objects.
   It will be changed later
*/

  array<typename gf_mesh<cyclic_lattice>::mesh_point_t, 1> r_arr(rmesh.size());
    auto r_iter = rmesh.begin();
    for (auto idx : range(0, rmesh.size())) {
        auto r  = *r_iter;
        r_arr(idx) = r;
        r_iter++;
    }
  auto _ = all_t{};
#pragma omp parallel for
  for(int idx_r = 0; idx_r < rmesh.size(); idx_r++){
    auto r = r_arr(idx_r);

    auto delta_t = make_gf<imtime>(tmesh, delta_tr_out.target());

    auto Gamma_pp_dyn_t = Gamma_pp_dyn_tr[_, r];
    auto F_t = F_tr[_, r];
    
    for (const auto t : tmesh) {
      for (auto [A, a, B, b] : Gamma_pp_dyn_tr.target_indices())
        delta_t[t](a, b) += -Gamma_pp_dyn_t[t](A, a, B, b) * F_t[t](A, B);
    }
    delta_tr_out[_, r] = delta_t;
  }
  t_dynamic_product.stop();
  
  // FIXME
  // This raises warnings when used with random delta input, e.g. eigenvalue finder
  t_fft_delta.start();
  //auto delta_wk_out = make_gf_from_fourier<0, 1>(delta_tr_out);
  auto delta_wr_out = g_wr_from_g_tr(delta_tr_out);
  auto delta_wk_out = g_wk_from_g_wr(delta_wr_out);
  t_fft_delta.stop();

  // Constant part
  t_constant_product.start();
  auto delta_k_out = eliashberg_constant_gamma_f_product(Gamma_pp_const_r, F_tr);
  t_constant_product.stop();

  // Combine dynamic and constant part
  t_combine.start();
  for (const auto [w , k]: delta_wk_out.mesh())
      delta_wk_out[w, k] += delta_k_out[k];
  t_combine.stop();

  std::cout << "g_delta_g_product:\t" << double(t_g_delta_g_product) << "\n" \
              << "fft_F:\t" << double(t_fft_F) << "\n" \
              << "dynamic_product:\t" << double(t_dynamic_product) << "\n" \
              << "fft_delta:\t" << double(t_fft_delta) << "\n" \
              << "constant_product:\t" << double(t_constant_product) << "\n" \
              << "combin:\t" << double(t_combine) << "\n" ;

  return delta_wk_out;
}

// optimized version if there is only a constant term
g_wk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r,
                                        g_wk_vt g_wk, g_wk_vt delta_wk) {

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);
  auto F_tr = make_gf_from_fourier<0, 1>(F_wk);

  auto delta_k_out = eliashberg_constant_gamma_f_product(Gamma_pp_const_r, F_tr);

  auto delta_wk_out = make_gf(F_wk);
  delta_wk_out *= 0.;

  for (const auto [w , k]: delta_wk_out.mesh())
      delta_wk_out[w, k] += delta_k_out[k];

  return delta_wk_out;
}

chi_wk_t gamma_PP_spin_charge(chi_wk_vt chi_c, chi_wk_vt chi_s, \
        array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s, \
        double charge_factor, double spin_factor) {

  using scalar_t = chi_wk_t::scalar_t;

  size_t nb = chi_c.target_shape()[0];

  auto Gamma_pp_wk = make_gf(chi_c);
  Gamma_pp_wk *= 0;

  // PH grouping of the vertex, from cc+cc+, permuting the last two indices.
  auto U_c_matrix = make_matrix_view(group_indices_view(U_c, {0, 1}, {3, 2}));
  auto U_s_matrix = make_matrix_view(group_indices_view(U_s, {0, 1}, {3, 2}));

  auto meshes_mpi = mpi_view(Gamma_pp_wk.mesh());

#pragma omp parallel for
  for (int idx = 0; idx < meshes_mpi.size(); idx++){
      auto &[w, k] = meshes_mpi(idx);

      array<scalar_t, 4> Gamma_pp_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};
      array<scalar_t, 4> chi_c_arr{chi_c[w, k], memory_layout_t<4>{0, 1, 2, 3}};
      array<scalar_t, 4> chi_s_arr{chi_s[w, k], memory_layout_t<4>{0, 1, 2, 3}};

      // PH grouping of the vertex, from cc+cc+, permuting the last two indices.
      auto Gamma_pp_matrix = make_matrix_view(group_indices_view(Gamma_pp_arr, {0, 1}, {3, 2}));
      // PH grouping of the susceptibilites, from c+cc+c, permuting the last two indices.
      auto chi_c_matrix = make_matrix_view(group_indices_view(chi_c_arr, {0, 1}, {3, 2}));
      auto chi_s_matrix = make_matrix_view(group_indices_view(chi_s_arr, {0, 1}, {3, 2}));

      Gamma_pp_matrix = charge_factor * U_c_matrix * chi_c_matrix * U_c_matrix \
                      + spin_factor * U_s_matrix * chi_s_matrix * U_s_matrix \
                      + 0.5 * (U_s_matrix + U_c_matrix);

      Gamma_pp_wk[w, k] = Gamma_pp_arr;
  }
  Gamma_pp_wk = mpi_all_reduce(Gamma_pp_wk);

  return Gamma_pp_wk;
}
  
chi_wk_t gamma_PP_singlet(chi_wk_vt chi_c, chi_wk_vt chi_s, \
        array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s) {

  auto Gamma_pp_wk = gamma_PP_spin_charge(chi_c, chi_s, U_c, U_s, -0.5, 1.5);
  return Gamma_pp_wk;
}

chi_wk_t gamma_PP_triplet(chi_wk_vt chi_c, chi_wk_vt chi_s, \
        array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s) {

  auto Gamma_pp_wk = gamma_PP_spin_charge(chi_c, chi_s, U_c, U_s, -0.5, -0.5);
  return Gamma_pp_wk;
}

} // namespace triqs_tprf
