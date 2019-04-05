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

namespace tprf {

// Helper function computing F = GG \Delta
  
gk_iw_t eliashberg_g_delta_g_product(gk_iw_vt g_wk, gk_iw_vt delta_wk) {

  auto [wmesh, kmesh] = delta_wk.mesh();
  auto gf_wmesh = std::get<0>(g_wk.mesh());

  if (wmesh.size() > gf_wmesh.size())
      TRIQS_RUNTIME_ERROR << "The size of the Matsubara frequency mesh of the Green's function"
          " (" << gf_wmesh.size() << ") must be atleast the size of the mesh of Delta (" <<
          wmesh.size() << ").";

  auto F_wk = make_gf(delta_wk);
  F_wk *= 0.;

  for (const auto [w, k] : delta_wk.mesh())
    for (auto [A, B] : F_wk.target_indices())
      for (auto [c, d] : delta_wk.target_indices())
        F_wk[w, k](A, B) +=
            g_wk[w, k](A, c) * g_wk[-w, -k](B, d) * delta_wk[w, k](c, d);

  return F_wk;
}

gk_iw_t eliashberg_product(chi_wk_vt Gamma_pp, gk_iw_vt g_wk,
                       gk_iw_vt delta_wk) {

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
      for (auto [A, a, b, B] : Gamma_pp.target_indices())
        delta_wk_out[w, k](a, b) +=
            Gamma_pp[w-n, k - q](A, a, b, B) * F_wk[n, q](A, B);

  delta_wk_out /= -(wmesh.domain().beta * kmesh.size());
  
  return delta_wk_out;
}

std::tuple<chi_wk_vt, chi_k_vt> split_into_dynamic_and_constant(chi_wk_vt Gamma_pp){

  auto _ = all_t{};
  auto [wmesh, kmesh] = Gamma_pp.mesh();
    
  // Fit infinite frequency value
  auto Gamma_pp_dyn = make_gf(Gamma_pp);

  chi_k_vt Gamma_pp_const_k = make_gf(kmesh, Gamma_pp.target());

  for (const auto k : kmesh) {
    auto Gamma_w = Gamma_pp[_, k];
    auto [tail, err] = fit_tail(Gamma_w);
    for (auto [a, b, c, d] : Gamma_pp.target_indices())
      Gamma_pp_const_k[k](a, b, c, d) = tail(0, a, b, c, d);
    for( const auto w : wmesh ) Gamma_pp_dyn[w, k] = Gamma_pp[w, k] - Gamma_pp_const_k[k];
  }

    return {Gamma_pp_dyn, Gamma_pp_const_k};
}


gk_iw_t eliashberg_product_fft(chi_wk_vt Gamma_pp,
                                gk_iw_vt g_wk, gk_iw_vt delta_wk) {

  auto _ = all_t{};

  auto F_wk = eliashberg_g_delta_g_product(g_wk, delta_wk);
  auto F_tr = make_gf_from_fourier<0, 1>(F_wk);

  auto [Gamma_pp_dyn, Gamma_pp_const_k] = split_into_dynamic_and_constant(Gamma_pp);

  auto Gamma_pp_tr = make_gf_from_fourier<0, 1>(Gamma_pp_dyn);

  auto delta_tr_out = make_gf(F_tr.mesh(), delta_wk.target());
  delta_tr_out *= 0.;

  for (const auto [t, r] : delta_tr_out.mesh()) {
    for (auto [A, a, b, B] : Gamma_pp.target_indices())
      delta_tr_out[t, r](a, b) += -Gamma_pp_tr(t, r)(A, a, b, B) * F_tr(t, r)(A, B);
  }
  
  auto delta_wk_out = make_gf_from_fourier<0, 1>(delta_tr_out);


  auto F_k = make_gf(Gamma_pp_const_k.mesh(), delta_wk.target());
  F_k *= 0.;

  for (const auto [w, k] : F_wk.mesh()) {
      for (auto [a, b] : delta_wk.target_indices())
          F_k[k](a, b) += F_wk[w, k](a, b);
  }

  auto F_r = make_gf_from_fourier<0>(F_k);
  auto Gamma_pp_const_r = make_gf_from_fourier<0>(Gamma_pp_const_k);
  
  auto delta_r_out = make_gf(F_r.mesh(), delta_wk.target());
  delta_r_out *= 0.;

  for (const auto r : delta_r_out.mesh()) {
    for (auto [A, a, b, B] : Gamma_pp.target_indices())
        delta_r_out[r](a, b) += -Gamma_pp_const_r(r)(A, a, b, B) * F_r(r)(A, B);
  }

  auto delta_k_out = make_gf_from_fourier<0>(delta_r_out);
  delta_k_out *= 1. / std::get<0>(delta_wk.mesh()).domain().beta;

  for (const auto [w , k]: delta_wk_out.mesh())
      delta_wk_out[w, k] += delta_k_out[k];

  return delta_wk_out;
}
  
chi_wk_t gamma_PP_singlet(chi_wk_vt chi_c, chi_wk_vt chi_s, \
        array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s) {

  auto [wmesh, kmesh] = chi_c.mesh();

  auto Gamma_pp_wk = make_gf(chi_c);
  Gamma_pp_wk *= 0;

  for (const auto [w, k] : Gamma_pp_wk.mesh())
    for (auto [a, b, c, d] : Gamma_pp_wk.target_indices())
      for (auto [A, B, C, D] : chi_c.target_indices())
        Gamma_pp_wk[w,k](a, b, c, d) += 
                1.5 * U_s(a, b, A, B) * chi_s[w, k](B, A, C, D) * U_s(D, C, c, d) \
                - 0.5 * U_c(a, b, A, B) * chi_c[w, k](B, A, C, D) * U_c(D, C, c, d) \
                + 0.5 * (U_s(a, b, c, d) + U_c(a, b, c, d));

  return Gamma_pp_wk;
}

chi_wk_t gamma_PP_triplet(chi_wk_vt chi_c, chi_wk_vt chi_s, \
        array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s) {

  auto [wmesh, kmesh] = chi_c.mesh();

  auto Gamma_pp_wk = make_gf(chi_c);
  Gamma_pp_wk *= 0;

  for (const auto [w, k] : Gamma_pp_wk.mesh())
    for (auto [a, b, c, d] : Gamma_pp_wk.target_indices())
      for (auto [A, B, C, D] : chi_c.target_indices())
        Gamma_pp_wk[w,k](a, b, c, d) += 
                - 0.5 * U_s(a, b, A, B) * chi_s[w, k](B, A, C, D) * U_s(D, C, c, d) \
                - 0.5 * U_c(a, b, A, B) * chi_c[w, k](B, A, C, D) * U_c(D, C, c, d) \
                + 0.5 * (U_s(a, b, c, d) + U_c(a, b, c, d));

  return Gamma_pp_wk;
}

} // namespace tprf
