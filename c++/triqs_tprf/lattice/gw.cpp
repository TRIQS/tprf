/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, The Simons Foundation
 * Authors: H. U.R. Strand, Y. in 't Veld, M. Rösner
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
#include "lattice_utility.hpp"
#include "../mpi.hpp"

// -- For parallell Fourier transform routines
#include "gf.hpp"
#include "chi_imtime.hpp"

namespace triqs_tprf {

  e_k_t rho_k_from_g_wk(g_wk_cvt g_wk) {

    auto _     = all_t{};
    auto kmesh = std::get<1>(g_wk.mesh());

    e_k_t rho_k(kmesh, g_wk.target_shape());
    rho_k() = 0.0;

    auto arr = mpi_view(kmesh);
#pragma omp parallel for
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &k = arr(idx);
      auto g_w  = g_wk(_, k);
      rho_k[k] = density(g_w);
    }
  
    rho_k = mpi::all_reduce(rho_k);
    return rho_k;
  }

  e_k_t rho_k_from_g_wk(g_Dwk_cvt g_wk) {

    auto _     = all_t{};
    auto wmesh = std::get<0>(g_wk.mesh());
    auto kmesh = std::get<1>(g_wk.mesh());

    e_k_t rho_k(kmesh, g_wk.target_shape());
    rho_k() = 0.0;

    auto arr = mpi_view(kmesh);
#pragma omp parallel for
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &k = arr(idx);
      auto g_w = make_gf<dlr_imfreq>({wmesh}, g_wk.target());
      g_w  = g_wk(_, k);
      rho_k[k] = density(dlr_coeffs_from_dlr_imfreq(g_w));
    }
  
    rho_k = mpi::all_reduce(rho_k);
    return rho_k;
  }
  
  template<typename W_t, typename g_t>
  auto gw_dynamic_sigma_impl(W_t W_tr, g_t g_tr) {

    auto Wtm = std::get<0>(W_tr.mesh());
    auto gtm = std::get<0>(g_tr.mesh());

    if (Wtm.size() != gtm.size() || Wtm.domain().beta != gtm.domain().beta) TRIQS_RUNTIME_ERROR << "gw_sigma_tr: tau meshes are not the same.\n";

    if (Wtm.domain().statistic != Boson || gtm.domain().statistic != Fermion) TRIQS_RUNTIME_ERROR << "gw_sigma_tr: statistics are incorrect.\n";

    if (std::get<1>(W_tr.mesh()) != std::get<1>(g_tr.mesh())) TRIQS_RUNTIME_ERROR << "gw_sigma_tr: real-space meshes are not the same.\n";

    auto sigma_tr = make_gf(g_tr);
    sigma_tr()    = 0.0;

    auto arr = mpi_view(g_tr.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[t, r] = arr(idx);

    for (const auto &[a, b, c, d] : W_tr.target_indices()) { sigma_tr[t, r](a, b) += -W_tr[t, r](a, c, d, b) * g_tr[t, r](c, d); }
  }

  sigma_tr = mpi::all_reduce(sigma_tr);
  return sigma_tr;
  }

  g_tr_t gw_dynamic_sigma(chi_tr_cvt W_tr, g_tr_cvt g_tr) {
    return gw_dynamic_sigma_impl(W_tr, g_tr);
  }

  g_Dtr_t gw_dynamic_sigma(chi_Dtr_cvt W_tr, g_Dtr_cvt g_tr) {
    return gw_dynamic_sigma_impl(W_tr, g_tr);
  }
  

  e_r_t hartree_sigma(chi_k_cvt v_k, e_r_cvt rho_r) {

  e_r_t sigma_r(rho_r.mesh(), rho_r.target_shape());
  sigma_r() = 0.0;

  for (auto const &[a, b, c, d] : v_k.target_indices()) {
    sigma_r[{0, 0, 0}](a, b) += v_k[{0, 0, 0}](a, b, c, d) * rho_r[{0, 0, 0}](c, d);
  }

  return sigma_r;
  }

  e_k_t hartree_sigma(chi_k_cvt v_k, g_wk_cvt g_wk) {

  if (v_k.mesh() != std::get<1>(g_wk.mesh())) TRIQS_RUNTIME_ERROR << "hartree_sigma: k-space meshes are not the same.\n";

  auto _     = all_t{};
  auto kmesh = std::get<1>(g_wk.mesh());

  e_k_t sigma_k(kmesh, g_wk.target_shape());
  sigma_k() = 0.0;

  auto arr = mpi_view(kmesh);
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &k = arr(idx);

    for (auto const &q : kmesh) {

      auto g_w  = g_wk(_, k + q);
      auto dens = density(g_w);

      for (auto const &[a, b, c, d] : v_k.target_indices()) { sigma_k[k](a, b) += v_k[q](a, b, c, d) * dens(c, d) / kmesh.size(); }
    }
  }
  sigma_k = mpi::all_reduce(sigma_k);
  return sigma_k;
  }

  e_r_t fock_sigma(chi_r_cvt v_r, e_r_cvt rho_r) {

  if (v_r.mesh() != rho_r.mesh()) TRIQS_RUNTIME_ERROR << "fock_sigma: r-space meshes are not the same.\n";

  auto rmesh = rho_r.mesh();

  e_r_t sigma_r(rmesh, rho_r.target_shape());
  sigma_r() = 0.0;

  auto arr = mpi_view(rmesh);
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &r = arr(idx);

    for (auto const &[a, b, c, d] : v_r.target_indices()) {
      sigma_r[r](a, b) += -v_r[r](a, c, d, b) * rho_r[r](d, c);
    }
  }
  sigma_r = mpi::all_reduce(sigma_r);
  return sigma_r;
  }

  e_k_t fock_sigma(chi_k_cvt v_k, g_wk_cvt g_wk) {

  if (v_k.mesh() != std::get<1>(g_wk.mesh())) TRIQS_RUNTIME_ERROR << "fock_sigma: k-space meshes are not the same.\n";

  auto _     = all_t{};
  auto kmesh = std::get<1>(g_wk.mesh());

  e_k_t sigma_k(kmesh, g_wk.target_shape());
  sigma_k() = 0.0;

  auto arr = mpi_view(kmesh);
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &k = arr(idx);

    for (auto const &q : kmesh) {

      auto g_w  = g_wk(_, k + q);
      auto dens = density(g_w);

      for (auto const &[a, b, c, d] : v_k.target_indices()) { sigma_k[k](a, b) += -v_k[q](a, c, d, b) * dens(d, c) / kmesh.size(); }
    }
  }
  sigma_k = mpi::all_reduce(sigma_k);
  return sigma_k;
  }

  e_k_t gw_sigma(chi_k_cvt v_k, g_wk_cvt g_wk) {
    return fock_sigma(v_k, g_wk);
  }

  template<typename W_t, typename g_t>
  auto gw_sigma_impl(W_t W_wk, g_t g_wk) {

  auto Wwm = std::get<0>(W_wk.mesh());
  auto gwm = std::get<0>(g_wk.mesh());

  if (Wwm.domain().beta != gwm.domain().beta)
    TRIQS_RUNTIME_ERROR << "gw_sigma: inverse temperatures are not the same.\n";
  if (Wwm.domain().statistic != Boson || gwm.domain().statistic != Fermion)
    TRIQS_RUNTIME_ERROR << "gw_sigma: statistics are incorrect.\n";
  if (std::get<1>(W_wk.mesh()) != std::get<1>(g_wk.mesh()))
    TRIQS_RUNTIME_ERROR << "gw_sigma: k-space meshes are not the same.\n";

  auto [W_dyn_wk, W_const_k] = split_into_dynamic_wk_and_constant_k(W_wk);

  // Dynamic GW self energy
  //auto g_tr = make_gf_from_fourier<0, 1>(g_wk); // Fixme! Use parallell transform
  auto g_wr = fourier_wk_to_wr(g_wk);
  auto g_tr = fourier_wr_to_tr(g_wr);
  
  //auto W_dyn_tr = make_gf_from_fourier<0, 1>(W_dyn_wk); // Fixme! Use parallell transform
  auto W_dyn_wr = chi_wr_from_chi_wk(W_dyn_wk);
  auto W_dyn_tr = chi_tr_from_chi_wr(W_dyn_wr);
  
  auto sigma_dyn_tr = gw_dynamic_sigma(W_dyn_tr, g_tr);

  // Static Fock part

  //auto sigma_fock_k = fock_sigma(W_const_k, g_wk); // Has N_k^2 scaling, not fast..

  auto W_const_r = make_gf_from_fourier(W_const_k);
  auto rho_k = rho_k_from_g_wk(g_wk);
  auto rho_r = make_gf_from_fourier(rho_k);
  auto sigma_fock_r = fock_sigma(W_const_r, rho_r);
  auto sigma_fock_k = make_gf_from_fourier(sigma_fock_r);

  // Add dynamic and static parts
  auto _ = all_t{};
  auto sigma_wk = make_gf_from_fourier<0, 1>(sigma_dyn_tr); // Fixme! Use parallell transform
  for (auto const &w : gwm) sigma_wk[w, _] += sigma_fock_k; // Single loop with no work per w, no need to parallellize over mpi

  return sigma_wk;
  }

  g_wk_t gw_sigma(chi_wk_cvt W_wk, g_wk_cvt g_wk) {
    return gw_sigma_impl(W_wk, g_wk);
  }

  /*
  g_Dwk_t gw_sigma(chi_Dwk_cvt W_wk, g_Dwk_cvt g_wk) {
    return gw_sigma_impl(W_wk, g_wk);
  }
  */
  
  // ----------------------------------------------------
  // g0w_sigma via spectral representation
  // dynamic part ...

  g_f_t g0w_dynamic_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, mesh::brzone::point_t kpoint) {

  if (std::get<1>(W_fk.mesh()) != e_k.mesh()) TRIQS_RUNTIME_ERROR << "g0w_sigma: k-space meshes are not the same.\n";
  if (e_k.mesh() != v_k.mesh()) TRIQS_RUNTIME_ERROR << "g0w_sigma: k-space meshes are not the same.\n";

  auto fmesh = std::get<0>(W_fk.mesh());
  auto kmesh = e_k.mesh();
  int nb     = e_k.target().shape()[0];

  std::complex<double> idelta(0.0, delta);

  g_f_t sigma_f(fmesh, e_k.target_shape());
  sigma_f() = 0.0;

  for (auto const &q : kmesh) {
    auto qpoint = mesh::brzone::point_t {q};
    auto kpqpoint = qpoint + kpoint;
    auto kpqvec = std::array<double, 3>{kpqpoint(0), kpqpoint(1), kpqpoint(2)};

    array<std::complex<double>, 2> e_kq_mat(e_k(kpqvec) - mu);
    auto eig_kq = linalg::eigenelements(e_kq_mat);
    auto ekq    = eig_kq.first;
    auto Ukq    = eig_kq.second;

    for (int l : range(nb)) {

      for (auto const &f : fmesh) {

        for (auto const &fp : fmesh) {

          auto num = bose(fp * beta) + fermi(ekq(l) * beta);
          auto den = f + idelta + fp - ekq(l);

          for (const auto &[a, b] : sigma_f.target_indices()) { 
            auto W_spec = -1.0 / M_PI * (W_fk[fp, q](a, a, b, b) - v_k[q](a, a, b, b)).imag();
            sigma_f[f](a, b) = sigma_f[f](a, b) 
                  + Ukq(a, l) * dagger(Ukq)(l, b) * (W_spec * num / den * fmesh.delta() / kmesh.size());            
          }
        }
      }
    }
  }
  
  return sigma_f;
  }

  g_fk_t g0w_dynamic_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, gf_mesh<brzone> kmesh) {

  if (std::get<1>(W_fk.mesh()) != e_k.mesh()) TRIQS_RUNTIME_ERROR << "g0w_sigma: k-space meshes are not the same.\n";
  if (e_k.mesh() != v_k.mesh()) TRIQS_RUNTIME_ERROR << "g0w_sigma: k-space meshes are not the same.\n";

  auto fmesh = std::get<0>(W_fk.mesh());

  g_fk_t sigma_fk({fmesh, kmesh}, e_k.target_shape());
  sigma_fk() = 0.0;

  auto arr = mpi_view(kmesh);
#pragma omp parallel for 
  for (unsigned int kidx = 0; kidx < arr.size(); kidx++) {
    auto &k = arr(kidx);
    auto kpoint = mesh::brzone::point_t {k};
    auto sigma_f = g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta, kpoint);
    for (auto const &f : fmesh) { sigma_fk[f,k] = sigma_f[f]; }
  }

  sigma_fk = mpi::all_reduce(sigma_fk);
  return sigma_fk;
  }

  g_fk_t g0w_dynamic_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta) {
  return g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta, e_k.mesh());
  }

  // static part ...

  array<std::complex<double>, 2> g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k, mesh::brzone::point_t kpoint) {
  
  if (e_k.mesh() != v_k.mesh()) TRIQS_RUNTIME_ERROR << "g0w_sigma: k-space meshes are not the same.\n";

  auto kmesh = e_k.mesh();
  int nb     = e_k.target().shape()[0];

  array<std::complex<double>, 2> sigma_k(nb ,nb);
  sigma_k() = 0.0;

  for (auto const &q : kmesh) {
    auto qpoint = mesh::brzone::point_t {q};
    auto kpqpoint = qpoint + kpoint;
    auto kpqvec = std::array<double, 3>{kpqpoint(0), kpqpoint(1), kpqpoint(2)};

    array<std::complex<double>, 2> e_kq_mat(e_k(kpqvec) - mu);
    auto eig_kq = linalg::eigenelements(e_kq_mat);
    auto ekq    = eig_kq.first;
    auto Ukq    = eig_kq.second;

    for (int l : range(nb)) {
      for (int a : range(nb)) {
        for (int b : range(nb)) {
          sigma_k(a, b) = sigma_k(a, b) - Ukq(a, l) * dagger(Ukq)(l, b) * v_k[q](a, a, b, b) * fermi(ekq(l) * beta) / kmesh.size();
        }
      }
    }
  }

  return sigma_k;
  }

  e_k_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k, gf_mesh<brzone> kmesh) {

  e_k_t sigma_k(kmesh, e_k.target_shape());
  sigma_k() = 0.0;

  auto arr = mpi_view(kmesh);
#pragma omp parallel for 
  for (unsigned int kidx = 0; kidx < arr.size(); kidx++) {
    auto &k = arr(kidx);
    auto kpoint = mesh::brzone::point_t {k};
    sigma_k[k] = g0w_sigma(mu, beta, e_k, v_k, kpoint);
  }

  sigma_k = mpi::all_reduce(sigma_k);
  return sigma_k;
  }

  e_k_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k) {
  return g0w_sigma(mu, beta, e_k, v_k, e_k.mesh());
  }

  // dynamic and static parts ...

  g_f_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, mesh::brzone::point_t kpoint) {
  auto sigma_stat_k = g0w_sigma(mu, beta, e_k, v_k, kpoint);
  auto sigma_dyn_fk = g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta, kpoint);

  auto fmesh = std::get<0>(W_fk.mesh());
  g_f_t sigma_fk(fmesh, e_k.target_shape());
  sigma_fk() = 0.0;
  for (auto const &f : fmesh) {
    sigma_fk[f] = sigma_stat_k + sigma_dyn_fk[f];
  }

  return sigma_fk;
  }

  g_fk_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, gf_mesh<brzone> kmesh) {
  auto sigma_stat_k = g0w_sigma(mu, beta, e_k, v_k, kmesh);
  auto sigma_dyn_fk = g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta, kmesh);
  auto sigma_fk = add_dynamic_fk_and_static_k(sigma_dyn_fk, sigma_stat_k);
  return sigma_fk;
  }
  
  g_fk_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta) {
  return g0w_sigma(mu, beta, e_k, W_fk, v_k, delta, e_k.mesh());
  }
} // namespace triqs_tprf
