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

chi_wk_t dynamical_screened_interaction_W_wk_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k) {

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
 


std::tuple<chi_wk_t, chi_k_t> split_into_dynamic_wk_and_constant_k(chi_wk_cvt W_wk) {

  auto _ = all_t{};
  auto wmesh = std::get<0>(W_wk.mesh());
  auto kmesh = std::get<1>(W_wk.mesh());
    
  chi_wk_t W_dyn_wk(W_wk.mesh(), W_wk.target_shape());
  chi_k_t W_const_k(kmesh, W_wk.target_shape());

  for (auto const &k : kmesh) {
    auto Gamma_w = W_wk[_, k];
    auto tail = std::get<0>(fit_tail(Gamma_w));

    for (auto [a, b, c, d] : W_wk.target_indices())
      W_const_k[k](a, b, c, d) = tail(0, a, b, c, d);

    for (auto const &w : wmesh) 
      W_dyn_wk[w, k] = W_wk[w, k] - W_const_k[k];
  }

  return {W_dyn_wk, W_const_k};
}


 
g_tr_t gw_sigma(chi_tr_cvt W_tr, g_tr_cvt g_tr) {

  auto Wtm = std::get<0>(W_tr.mesh());
  auto gtm = std::get<0>(g_tr.mesh());
  
  if( Wtm.size() != gtm.size() || Wtm.domain().beta != gtm.domain().beta )
    TRIQS_RUNTIME_ERROR << "gw_sigma_tr: tau meshes are not the same.\n";

  if( Wtm.domain().statistic != Boson || gtm.domain().statistic != Fermion )
    TRIQS_RUNTIME_ERROR << "gw_sigma_tr: statistics are incorrect.\n";
  
  if( std::get<1>(W_tr.mesh()) != std::get<1>(g_tr.mesh()) )
    TRIQS_RUNTIME_ERROR << "gw_sigma_tr: real-space meshes are not the same.\n";
  
  auto sigma_tr = make_gf(g_tr);
  sigma_tr *= 0.;

  auto arr = mpi_view(g_tr.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[t, r] = arr(idx);

    for (const auto &[a, b, c, d] : W_tr.target_indices()) {
      sigma_tr[t, r](a, b) += W_tr[t, r](a, b, c, d) * g_tr[t, r](c, d);
    }
  }

  sigma_tr = mpi::all_reduce(sigma_tr);
  return sigma_tr;
}

e_k_t gw_sigma(chi_k_cvt v_k, g_wk_cvt g_wk){

  if( v_k.mesh() != std::get<1>(g_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "gw_sigma: k-space meshes are not the same.\n";

  auto _ = all_t{};
  auto kmesh = std::get<1>(g_wk.mesh());

  e_k_t sigma_k(kmesh, g_wk.target_shape());
  
  auto arr = mpi_view(kmesh);
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto k = arr(idx);

    for (auto const &q : kmesh) {
      
      auto density = triqs::gfs::density(g_wk[_,k+q]);

      for (const auto &[a, b] : sigma_k.target_indices()) {
        sigma_k[k](a, b) += - v_k[q](a, b, a, b) * density(a,b) / kmesh.size();
      }
    }
  }
  sigma_k = mpi::all_reduce(sigma_k);
  return sigma_k;
}

g_wk_t gw_sigma(chi_wk_cvt W_wk, g_wk_cvt g_wk) {

  auto Wwm = std::get<0>(W_wk.mesh());
  auto gwm = std::get<0>(g_wk.mesh());
  
  if( Wwm.domain().beta != gwm.domain().beta )
    TRIQS_RUNTIME_ERROR << "gw_sigma: inverse temperatures are not the same.\n";

  if( Wwm.domain().statistic != Boson || gwm.domain().statistic != Fermion )
    TRIQS_RUNTIME_ERROR << "gw_sigma: statistics are incorrect.\n";

  if( std::get<1>(W_wk.mesh()) != std::get<1>(g_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "gw_sigma: k-space meshes are not the same.\n";


  auto [W_dyn_wk, W_const_k] = split_into_dynamic_wk_and_constant_k(W_wk);

  //Dynamic part
  auto g_tr = make_gf_from_fourier<0, 1>(g_wk);
  auto W_tr = make_gf_from_fourier<0, 1>(W_dyn_wk);
  auto sigma_tr = gw_sigma(W_tr, g_tr);
  auto sigma_dyn_wk = make_gf_from_fourier<0, 1>(sigma_tr);

  //Static part
  auto sigma_stat_k = gw_sigma(W_const_k, g_wk);

  //Add dynamic and static parts
  g_wk_t sigma_wk(g_wk.mesh(), g_wk.target_shape());

  auto arr = mpi_view(sigma_wk.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    for (const auto &[a, b] : sigma_wk.target_indices()) {
      sigma_wk[w, k](a, b) += sigma_dyn_wk[w, k](a, b) + sigma_stat_k[k](a, b);
    }
  }
  sigma_wk = mpi::all_reduce(sigma_wk);
  return sigma_wk;
}

// ----------------------------------------------------
// helper functions

double fermi2(double e) { return 1. / (exp(e) + 1.); }
double bose2(double e)  { return 1. / (exp(e) - 1.); }

// ----------------------------------------------------
// gw_sigma_fk_g0w0_spectral

g_fk_t g0w_sigma(double mu, double beta, e_k_cvt e_k, 
                 gf_mesh<refreq> mesh, chi_fk_cvt W_fk, 
                 chi_k_cvt v_k, double delta) {

  auto fmesh = mesh;
  auto kmesh = e_k.mesh();
  int nb = e_k.target().shape()[0];
  
  std::complex<double> idelta(0.0, delta);
  
  g_fk_t sigma_fk({fmesh, kmesh}, e_k.target_shape());
  
  for (auto const & [ f, k ] : sigma_fk.mesh())
    sigma_fk[f, k] = 0.;
 
  g_fk_t WSpec_fk({fmesh, kmesh}, e_k.target_shape());
 
  for (auto const & [ f, k ] : WSpec_fk.mesh()) {
    for (int i : range(nb)) {
      for (int j : range(nb)) {
        WSpec_fk[f, k](i, j) = -1.0/3.141592653589793 * (W_fk[f, k](i, j, i, j) - v_k[k](i, j, i, j)).imag();
      }
    }
  }

  auto arr = mpi_view(kmesh);  
  #pragma omp parallel for shared(sigma_fk)
  for (unsigned int kidx = 0; kidx < arr.size(); kidx++) {
    auto k = arr(kidx);

    for (auto const &q : kmesh) {
    
      array<std::complex<double>, 2> e_kq_mat(e_k(k + q) - mu);
      auto eig_kq = linalg::eigenelements(e_kq_mat);
      auto ekq = eig_kq.first;
      auto Ukq = eig_kq.second;

      for (int l : range(nb)) {

        for (auto const &f : fmesh) {
            
              for (auto const &fp : fmesh) {
              
                auto num   = bose2(fp * beta) + fermi2(ekq(l) * beta);
                auto den   = f + idelta + fp - ekq(l);

                sigma_fk[f,k](a,b)
                  << sigma_fk[f,k](a,b) + Ukq(l, a) * dagger(Ukq)(b, l) * \
                      ( WSpec_fk[fp, q](a, b) * num / den * fmesh.delta() / kmesh.size() );

              }
                
              sigma_fk[f,k](a,b)
                << sigma_fk[f,k](a,b) - Ukq(l, a) * dagger(Ukq)(b, l) * \
                                        v_k[q](a, b, a, b) * fermi2(ekq(l) * beta) / kmesh.size();
            
         }
      }
    }
  } 

  sigma_fk = mpi::all_reduce(sigma_fk);
  return sigma_fk;
}

// ----------------------------------------------------
// gw_sigma_k_g0w0

e_k_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k) {
  auto kmesh = e_k.mesh();
  int nb = e_k.target().shape()[0];

  e_k_t sigma_k(kmesh, e_k.target_shape());

  for (auto const &k : sigma_k.mesh())
    sigma_k[k] = 0.; 
 
  auto arr = mpi_view(kmesh);
  #pragma omp parallel for
  for (unsigned int kidx = 0; kidx < arr.size(); kidx++) {
    auto k = arr(kidx);
    
    for (auto const &q : kmesh) {
      array<std::complex<double>, 2> e_kq_mat(e_k(k + q) - mu);
      auto eig_kq = linalg::eigenelements(e_kq_mat);
      auto ekq = eig_kq.first;
      auto Ukq = eig_kq.second;

      for (auto [a,b] : sigma_k.target_indices()){
        for (int l : range(nb)) {
          sigma_k[k](a,b) += - Ukq(l, a) * dagger(Ukq)(b, l) * \
                               v_k[q](a, b, a, b) * fermi2(ekq(l) * beta) / kmesh.size();
        }
      }
    }
  }

  sigma_k = mpi::all_reduce(sigma_k);
  return sigma_k;
}


} // namespace triqs_tprf
