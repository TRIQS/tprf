/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation and S. Käser
 * Authors: S. Käser, H. U.R. Strand
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

#include "eliashberg.hpp"
#include <omp.h>
#include "../mpi.hpp"

#include "gf.hpp"
#include "fourier.hpp"

namespace triqs_tprf {

// Helper function computing F = GG \Delta

g_wk_t eliashberg_g_delta_g_product(g_wk_vt g_wk, g_wk_vt delta_wk) {

  // Get rid of structured binding declarations in this file due to issue #11
  //auto [wmesh, kmesh] = delta_wk.mesh();
  auto wmesh = std::get<0>(delta_wk.mesh());
  auto kmesh = std::get<1>(delta_wk.mesh());

  auto wmesh_gf = std::get<0>(g_wk.mesh());

  if (wmesh.size() > wmesh_gf.size())
      TRIQS_RUNTIME_ERROR << "The size of the Matsubara frequency mesh of the Green's function"
          " (" << wmesh_gf.size() << ") must be atleast the size of the mesh of Delta (" <<
          wmesh.size() << ").";

  auto F_wk = make_gf(delta_wk);
  F_wk *= 0.;

  auto meshes_mpi = mpi_view(delta_wk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < meshes_mpi.size(); idx++){
      auto &[w, k] = meshes_mpi[idx];

    for (auto [d, c] : F_wk.target_indices()){
      for (auto [e, f] : delta_wk.target_indices()){
        F_wk[w, k](d, c) +=
        g_wk[w, k](c, f) * g_wk[-w, -k](d, e) * delta_wk[w, k](e, f);
      }
    }
  }

  F_wk = mpi::all_reduce(F_wk);

  return F_wk;
}

g_Dwk_t eliashberg_g_delta_g_product(g_Dwk_vt g_wk, g_Dwk_vt delta_wk) {

  // Get rid of structured binding declarations in this file due to issue #11
  //auto [wmesh, kmesh] = delta_wk.mesh();
  auto wmesh = std::get<0>(delta_wk.mesh());
  auto kmesh = std::get<1>(delta_wk.mesh());

  auto wmesh_gf = std::get<0>(g_wk.mesh());

  //auto cmesh = triqs::mesh::dlr_coeffs(wmesh);

  if (wmesh.size() > wmesh_gf.size())
      TRIQS_RUNTIME_ERROR << "The size of the Matsubara frequency mesh of the Green's function"
          " (" << wmesh_gf.size() << ") must be atleast the size of the mesh of Delta (" <<
          wmesh.size() << ").";

  auto F_wk = make_gf(delta_wk);
  F_wk *= 0.;

  auto _ = all_t{};

  auto meshes_mpi = mpi_view(delta_wk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < meshes_mpi.size(); idx++){
    auto &[w, k] = meshes_mpi(idx);

    // This does not work due to problem in cppdlr
    //auto g_Dck = dlr_coeffs_from_dlr_imfreq(g_wk[_,k]);
    //auto g_Dcmk = dlr_coeffs_from_dlr_imfreq(g_wk[_,-k]);
    //auto delta_Dck = dlr_coeffs_from_dlr_imfreq(delta_wk[_,k]);

    g_Dw_t g_w({wmesh}, g_wk.target_shape());
    g_w() = g_wk[_,k];
    auto g_Dck = dlr_coeffs_from_dlr_imfreq(g_w);
    g_Dw_t g_w2({wmesh}, g_wk.target_shape());
    g_w2() = g_wk[_,-k];
    auto g_Dcmk = dlr_coeffs_from_dlr_imfreq(g_w2);
    g_Dw_t delta_w({wmesh}, delta_wk.target_shape());
    delta_w() = delta_wk[_,k];
    auto delta_Dck = dlr_coeffs_from_dlr_imfreq(delta_w);

    for (auto [d, c] : F_wk.target_indices()) {
      for (auto [e, f] : delta_wk.target_indices()) {
        F_wk[w, k](d, c) += g_Dck(w)(c, f) * g_Dcmk(-w)(d, e) * delta_Dck(w)(e, f);
      }
    }
  }

  F_wk = mpi::all_reduce(F_wk);

  return F_wk;
}


g_wk_t eliashberg_product(chi_wk_vt Gamma_pp, g_wk_vt g_wk,
                       g_wk_vt delta_wk) {

  //auto [wmesh, kmesh] = delta_wk.mesh();
  auto wmesh = std::get<0>(delta_wk.mesh());
  auto kmesh = std::get<1>(delta_wk.mesh());

  auto gamma_wmesh = std::get<0>(Gamma_pp.mesh());

  if (2*wmesh.size() > gamma_wmesh.size())
      TRIQS_RUNTIME_ERROR << "The size of the Matsubara frequency mesh of Gamma"
          " (" << gamma_wmesh.size() << ") must be atleast TWICE the size of the mesh of Delta (" <<
          wmesh.size() << ").";

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);

  auto delta_wk_out = make_gf(delta_wk);
  delta_wk_out *= 0.;

  auto arr = mpi_view(kmesh);
#pragma omp parallel for
  for (int kidx = 0; kidx < arr.size(); kidx++) {
      auto &k = arr[kidx];
      for (auto w : wmesh) {
        for (auto [n, q] : delta_wk.mesh())
          for (auto [c, a, d, b] : Gamma_pp.target_indices())
            delta_wk_out[w, k](a, b) += -0.5 * Gamma_pp(w - n, k - q)(c, a, d, b) * F_wk[n, q](d, c);
      }
  }

  delta_wk_out /= (wmesh.beta() * kmesh.size());

  delta_wk_out = mpi::all_reduce(delta_wk_out);

  return delta_wk_out;
}

std::tuple<chi_tr_t, chi_r_t> dynamic_and_constant_to_tr(chi_wk_vt Gamma_pp_dyn_wk, chi_k_vt Gamma_pp_const_k) {

    auto Gamma_pp_dyn_wr = fourier_wk_to_wr_general_target(Gamma_pp_dyn_wk);
    auto Gamma_pp_dyn_tr = fourier_wr_to_tr_general_target(Gamma_pp_dyn_wr);

    auto Gamma_pp_const_r = make_gf_from_fourier<0>(Gamma_pp_const_k);

    return {Gamma_pp_dyn_tr, Gamma_pp_const_r}; 
}

std::tuple<chi_Dtr_t, chi_r_t> dynamic_and_constant_to_tr(chi_Dwk_vt Gamma_pp_dyn_wk, chi_k_vt Gamma_pp_const_k) {

    auto Gamma_pp_dyn_wr = fourier_wk_to_wr_general_target(Gamma_pp_dyn_wk);
    auto Gamma_pp_dyn_tr = fourier_Dwr_to_Dtr_general_target(Gamma_pp_dyn_wr);

    auto Gamma_pp_const_r = make_gf_from_fourier<0>(Gamma_pp_const_k);

    return {Gamma_pp_dyn_tr, Gamma_pp_const_r}; 
}

e_r_t eliashberg_constant_gamma_f_product(chi_r_vt Gamma_pp_const_r, g_tr_t F_tr) {

  auto _ = all_t{};

  auto delta_r_out = make_gf(std::get<1>(F_tr.mesh()), F_tr.target());
  delta_r_out *= 0.;

  for (auto r : std::get<1>(F_tr.mesh())) {
      auto F_t = F_tr[_, r];
      for (auto [c, a, d, b] : Gamma_pp_const_r.target_indices()) delta_r_out[r](a, b) += -0.5 * Gamma_pp_const_r[r](c, a, d, b) * F_t(0)(d, c);
  }

  return delta_r_out;
}

e_r_t eliashberg_constant_gamma_f_product(chi_r_vt Gamma_pp_const_r, g_Dtr_t F_tr) {

  auto _ = all_t{};
  auto tmesh = std::get<0>(F_tr.mesh());

  auto delta_r_out = make_gf(std::get<1>(F_tr.mesh()), F_tr.target());
  delta_r_out *= 0.;

  for (auto const &r : std::get<1>(F_tr.mesh())) {
    g_Dt_t F_t({tmesh}, F_tr.target_shape());
    F_t() = F_tr[_, r];
    auto F_Dc = dlr_coeffs_from_dlr_imtime(F_t);
    for (auto [c, a, d, b] : Gamma_pp_const_r.target_indices())
        delta_r_out[r](a, b) += -0.5 * Gamma_pp_const_r[r](c, a, d, b) * F_Dc(0)(d, c);
  }

  return delta_r_out;
}

template<typename delta_t, typename chi_t, typename F_t>  
delta_t eliashberg_dynamic_gamma_f_product_template(chi_t Gamma_pp_dyn_tr, F_t F_tr) {

  //auto [tmesh, rmesh] = F_tr.mesh();
  auto tmesh = std::get<0>(F_tr.mesh());
  auto rmesh = std::get<1>(F_tr.mesh());

  auto delta_tr_out = make_gf(F_tr);
  delta_tr_out *= 0.;

  auto tmesh_gamma = std::get<0>(Gamma_pp_dyn_tr.mesh());

  // Test if the tau meshs of delta and gamma are compatible. If not raise an error, because
  // it would lead to wrong results.
  if (tmesh.size() != tmesh_gamma.size()) 
      TRIQS_RUNTIME_ERROR << "The size of the imaginary time mesh of Gamma"
          " (" << tmesh_gamma.size() << ") must be the size of the mesh of Delta (" <<
          tmesh.size() << ").";

  auto meshes_mpi = mpi_view(F_tr.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < meshes_mpi.size(); idx++){
      auto &[t, r] = meshes_mpi[idx];

      for (auto [c, a, d, b] : Gamma_pp_dyn_tr.target_indices())
        delta_tr_out[t, r](a, b) += -0.5 * Gamma_pp_dyn_tr[t, r](c, a, d, b) * F_tr[t, r](d, c);
  }

  delta_tr_out = mpi::all_reduce(delta_tr_out);

  return delta_tr_out;
}

g_tr_t eliashberg_dynamic_gamma_f_product(chi_tr_vt Gamma_pp_dyn_tr, g_tr_vt F_tr) {
  return eliashberg_dynamic_gamma_f_product_template<g_tr_t, chi_tr_vt, g_tr_vt>(Gamma_pp_dyn_tr, F_tr);
}

g_Dtr_t eliashberg_dynamic_gamma_f_product(chi_Dtr_vt Gamma_pp_dyn_tr, g_Dtr_vt F_tr) {
  return eliashberg_dynamic_gamma_f_product_template<g_Dtr_t, chi_Dtr_vt, g_Dtr_vt>(Gamma_pp_dyn_tr, F_tr);
}


template<typename delta_out_t, typename chi_t, typename g_t>  
delta_out_t eliashberg_product_fft_template(chi_t Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r,
                                   g_t g_wk, g_t delta_wk) {

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);
  auto F_wr = fourier_wk_to_wr(F_wk);
  auto F_tr = fourier_wr_to_tr(F_wr);

  auto delta_tr_out = eliashberg_dynamic_gamma_f_product(Gamma_pp_dyn_tr, F_tr);
  auto delta_r_out = eliashberg_constant_gamma_f_product(Gamma_pp_const_r, F_tr);
  
  // FIXME
  // This raises warnings when used with random delta input, e.g. eigenvalue finder
  auto delta_wr_out = fourier_tr_to_wr(delta_tr_out);
  // Combine dynamic and constant part
  auto _ = all_t{};
  for (auto w : std::get<0>(delta_wr_out.mesh())) delta_wr_out[w, _] += delta_r_out;

  auto delta_wk_out = fourier_wr_to_wk(delta_wr_out);

  return delta_wk_out;
}

g_wk_t eliashberg_product_fft(chi_tr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt delta_wk) {
  return eliashberg_product_fft_template<g_wk_t, chi_tr_vt, g_wk_vt>(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk);
}

g_Dwk_t eliashberg_product_fft(chi_Dtr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_Dwk_vt g_wk, g_Dwk_vt delta_wk) {
  return eliashberg_product_fft_template<g_Dwk_t, chi_Dtr_vt, g_Dwk_vt>(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk);
}

// optimized version if there is only a constant term

template<typename delta_out_t, typename g_t>  
delta_out_t eliashberg_product_fft_constant_template(chi_r_vt Gamma_pp_const_r,
                                        g_t g_wk, g_t delta_wk) {

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);
  auto F_wr = fourier_wk_to_wr(F_wk);
  auto F_tr = fourier_wr_to_tr(F_wr);

  auto delta_r_out = eliashberg_constant_gamma_f_product(Gamma_pp_const_r, F_tr);
  auto delta_k_out = make_gf_from_fourier<0>(delta_r_out);

  auto delta_wk_out = make_gf(F_wk);
  delta_wk_out *= 0.;

  auto _ = all_t{};
  for (auto w : std::get<0>(delta_wk_out.mesh())) delta_wk_out[w, _] += delta_k_out;

  return delta_wk_out;
}

g_wk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt delta_wk) {
  return eliashberg_product_fft_constant_template<g_wk_t, g_wk_vt>(Gamma_pp_const_r, g_wk, delta_wk);
}

g_Dwk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_Dwk_vt g_wk, g_Dwk_vt delta_wk) {
  return eliashberg_product_fft_constant_template<g_Dwk_t, g_Dwk_vt>(Gamma_pp_const_r, g_wk, delta_wk);
}


chi_wk_t construct_phi_wk(chi_wk_vt chi, array_contiguous_view<std::complex<double>, 4> U) {

  using scalar_t = chi_wk_t::scalar_t;

  size_t nb = chi.target_shape()[0];

  auto phi_wk = make_gf(chi);
  phi_wk *= 0;

  // PH grouping of the vertex, from cc+cc+, permuting the last two indices.
  auto U_matrix = make_matrix_view(group_indices_view(U, idx_group<0, 1>, idx_group<3, 2>));

  auto meshes_mpi = mpi_view(phi_wk.mesh());

#pragma omp parallel for
  for (unsigned int idx = 0; idx < meshes_mpi.size(); idx++){
      auto &[w, k] = meshes_mpi[idx];

      array<scalar_t, 4> phi_arr{nb, nb, nb, nb};
      array<scalar_t, 4> chi_arr{chi[w, k]};

      // PH grouping of the vertex, from cc+cc+, permuting the last two indices.
      auto phi_matrix = make_matrix_view(group_indices_view(phi_arr, idx_group<0, 1>, idx_group<3, 2>));
      // PH grouping of the susceptibilites, from c+cc+c, permuting the last two indices.
      auto chi_matrix = make_matrix_view(group_indices_view(chi_arr, idx_group<0, 1>, idx_group<3, 2>));

      phi_matrix = U_matrix * chi_matrix * U_matrix;

      phi_wk[w, k] = phi_arr;
  }
  phi_wk = mpi::all_reduce(phi_wk);

  return phi_wk;
}

} // namespace triqs_tprf
