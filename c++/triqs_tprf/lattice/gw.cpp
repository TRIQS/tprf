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

#include <triqs/arrays/linalg/eigenelements.hpp>

#include "gw.hpp"
#include "common.hpp"
#include "../mpi.hpp"

namespace triqs_tprf {

chi_wk_t dynamical_screened_interaction_W_wk(chi_wk_cvt PI_wk, chi_k_cvt V_k) {

  if( std::get<1>(PI_wk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "retarded_screened_interaction_Wr_wk: k-space meshes are not the same\n";
  
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

    W = V * inverse(I - PI * V) - V;

    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

chi_wk_t dynamical_screened_interaction_W_wk_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k) {

  if( std::get<1>(chi_wk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "retarded_screened_interaction_Wr_wk: k-space meshes are not the same\n";
  
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

    W = V * chi * V;
    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}
  
g_tr_t gw_sigma_tr(chi_tr_cvt Wr_tr, g_tr_cvt g_tr) {

  auto Wtm = std::get<0>(Wr_tr.mesh());
  auto gtm = std::get<0>(g_tr.mesh());
  
  if( Wtm.size() != gtm.size() || Wtm.domain().beta != gtm.domain().beta )
    TRIQS_RUNTIME_ERROR << "gw_sigma_tr: tau meshes are not the same.\n";

  if( Wtm.domain().statistic != Boson || gtm.domain().statistic != Fermion )
    TRIQS_RUNTIME_ERROR << "gw_sigma_tr: statistics are incorrect.\n";
  
  if( std::get<1>(Wr_tr.mesh()) != std::get<1>(g_tr.mesh()) )
    TRIQS_RUNTIME_ERROR << "gw_sigma_tr: real-space meshes are not the same.\n";
  
  auto sigma_tr = make_gf(g_tr);
  sigma_tr *= 0.;

  // MPI and openMP parallell loop
  auto arr = mpi_view(g_tr.mesh());
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[t, r] = arr(idx);

    //for (const auto &[t, r] : g_tr.mesh()) {

    for (const auto &[a, b, c, d] : Wr_tr.target_indices()) {
      sigma_tr[t, r](a, b) += Wr_tr[t, r](a, b, c, d) * g_tr[t, r](c, d);
    }
  }

  sigma_tr = mpi::all_reduce(sigma_tr);
  return sigma_tr;
}

g_wk_t gw_sigma_wk_serial_fft(chi_wk_cvt Wr_wk, g_wk_cvt g_wk) {

  auto Wwm = std::get<0>(Wr_wk.mesh());
  auto gwm = std::get<0>(g_wk.mesh());
  
  if( Wwm.domain().beta != gwm.domain().beta )
    TRIQS_RUNTIME_ERROR << "gw_self_energy: inverse temperatures are not the same.\n";

  if( Wwm.domain().statistic != Boson || gwm.domain().statistic != Fermion )
    TRIQS_RUNTIME_ERROR << "gw_self_energy: statistics are incorrect.\n";

  if( std::get<1>(Wr_wk.mesh()) != std::get<1>(g_wk.mesh()) )
    TRIQS_RUNTIME_ERROR << "gw_self_energy: k-space meshes are not the same.\n";
  
  // TODO: parallellize fourier transforms
  auto g_tr = make_gf_from_fourier<0, 1>(g_wk);
  auto Wr_tr = make_gf_from_fourier<0, 1>(Wr_wk);
  auto sigma_tr = gw_sigma_tr(Wr_tr, g_tr);
  auto sigma_wk = make_gf_from_fourier<0, 1>(sigma_tr);
  return sigma_wk;
}

// ----------------------------------------------------
// helper functions

double fermi2(double e) { return 1. / (exp(e) + 1.); }
double bose2(double e)  { return 1. / (exp(e) - 1.); }

// ----------------------------------------------------
// gw_sigma_fk_g0w0_spectral

g_fk_t gw_sigma_fk_g0w0_spectral(double mu, double beta, e_k_cvt e_k, 
                                 gf_mesh<refreq> fmesh, chi_fk_cvt Wr_fk, 
                                 chi_k_cvt v_k, double delta) {

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
        WSpec_fk[f, k](i, j) = -1.0/3.141592653589793 * (Wr_fk[f, k](i, j, i, j) - v_k[k](i, j, i, j)).imag();
      }
    }
  }
  
  #pragma omp parallel for shared(sigma_fk)
  for (int kidx = 0; kidx < kmesh.size(); kidx++) {
  
    auto k_iter = kmesh.begin();
    k_iter += kidx;
    auto k = *k_iter;

    for (auto const &q : kmesh) {
    
      matrix<std::complex<double>> e_kq_mat(e_k(k + q) - mu);
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

  return sigma_fk;
}

} // namespace triqs_tprf
