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
#include "../mpi.hpp"

namespace triqs_tprf {

// ----------------------------------------------------
// dynamical_screened_interaction_W_fk

chi_fk_t dynamical_screened_interaction_W_fk(chi_fk_cvt pi0_fk, chi_k_cvt v_k) {

  if( std::get<1>(pi0_fk.mesh()) != v_k.mesh() )
    TRIQS_RUNTIME_ERROR << "dynamical_screened_interaction_W_fk: k-space meshes are not the same\n";
  
  auto W_fk = make_gf(pi0_fk);
  W_fk *= 0.;
  size_t nb = pi0_fk.target_shape()[0];

  using scalar_t = chi_fk_t::scalar_t;
  auto I = make_unit_matrix<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_fk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[f, k] = arr(idx);

    array<scalar_t, 4> v_arr{v_k[k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> pi0_arr{pi0_fk[f, k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};

    auto v   = make_matrix_view(group_indices_view(v_arr,   {0, 1}, {3, 2}));
    auto pi0 = make_matrix_view(group_indices_view(pi0_arr, {0, 1}, {3, 2}));
    auto W   = make_matrix_view(group_indices_view(W_arr,   {0, 1}, {3, 2}));

    W = v * inverse(I - pi0 * v);

    W_fk[f, k] = W_arr;

  }

  W_fk = mpi::all_reduce(W_fk);

  return W_fk;
  
}

// ----------------------------------------------------
// dynamical_screened_interaction_W_wk

chi_wk_t dynamical_screened_interaction_W_wk(chi_wk_cvt PI_wk, chi_k_cvt V_k) {

  if( std::get<1>(PI_wk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "retarded_screened_interaction_Wr_wk: k-space meshes are not the same\n";
  
  auto W_wk = make_gf(PI_wk);
  W_wk *= 0.;
  size_t nb = PI_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = make_unit_matrix<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_k[k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> PI_arr{PI_wk[w, k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};

    auto V = make_matrix_view(group_indices_view(V_arr, {0, 1}, {3, 2}));
    auto PI = make_matrix_view(group_indices_view(PI_arr, {0, 1}, {3, 2}));
    auto W = make_matrix_view(group_indices_view(W_arr, {0, 1}, {3, 2}));

    W = V * inverse(I - PI * V) - V;

    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

// ----------------------------------------------------
// dynamical_screened_interaction_W_wk_from_generalized_susceptibility

chi_wk_t dynamical_screened_interaction_W_wk_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k) {

  if( std::get<1>(chi_wk.mesh()) != V_k.mesh() )
    TRIQS_RUNTIME_ERROR << "retarded_screened_interaction_Wr_wk: k-space meshes are not the same\n";
  
  auto W_wk = make_gf(chi_wk);
  W_wk *= 0.;
  size_t nb = W_wk.target_shape()[0];

  using scalar_t = chi_wk_t::scalar_t;
  auto I = make_unit_matrix<scalar_t>(nb * nb);

  // MPI and openMP parallell loop
  auto arr = mpi_view(W_wk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 4> V_arr{V_k[k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> chi_arr{chi_wk[w, k], memory_layout_t<4>{0, 1, 2, 3}};
    array<scalar_t, 4> W_arr{nb, nb, nb, nb, memory_layout_t<4>{0, 1, 2, 3}};

    auto V = make_matrix_view(group_indices_view(V_arr, {0, 1}, {3, 2}));
    auto chi = make_matrix_view(group_indices_view(chi_arr, {0, 1}, {3, 2}));
    auto W = make_matrix_view(group_indices_view(W_arr, {0, 1}, {3, 2}));

    W = V * chi * V;
    W_wk[w, k] = W_arr;
  }

  W_wk = mpi::all_reduce(W_wk);
  return W_wk;
}

// ----------------------------------------------------
// helper functions

double fermi2(double e) { return 1. / (exp(e) + 1.); }
double bose2(double e)  { return 1. / (exp(e) - 1.); }

// ----------------------------------------------------
// gw_sigma_fk_g0w0_spectral

g_fk_t gw_sigma_fk_g0w0_spectral(double mu, double beta, h_k_cvt h_k, 
                                 gf_mesh<refreq> fmesh, chi_fk_cvt Wr_fk, 
                                 chi_k_cvt v_k, double delta) {

  auto kmesh = h_k.mesh();
  int nb = h_k.target().shape()[0];
  
  std::complex<double> idelta(0.0, delta);
  
  g_fk_t sigma_fk_dyn({fmesh, kmesh}, h_k.target_shape());
  g_fk_t sigma_fk_sta({fmesh, kmesh}, h_k.target_shape());
  g_fk_t sigma_fk({fmesh, kmesh}, h_k.target_shape());
  
//#ifdef TPRF_OMP
//  std::cout << "OMP ACTIVATED";
//#endif

  #pragma omp parallel for shared(sigma_fk_dyn, sigma_fk_sta)
  for (int kidx = 0; kidx < kmesh.size(); kidx++) {
  
    auto k_iter = kmesh.begin();
    k_iter += kidx;
    auto k = *k_iter;

    for (auto const &q : kmesh) {
    
      matrix<std::complex<double>> h_kq_mat(h_k(k + q) - mu);
      auto eig_kq = linalg::eigenelements(h_kq_mat);
      auto ekq = eig_kq.first;
      auto Ukq = eig_kq.second;
    
      for (int i : range(nb)) {
        for (int j : range(nb)) {
      
          for (auto const &f : fmesh) {

            for (int l : range(nb)) {
            
              sigma_fk_sta[f,k](j,i) -= fermi2(ekq(l) * beta) * v_k[q](i, j, i, j) * Ukq(l, i) * dagger(Ukq)(j, l);
      
              for (auto const &fp : fmesh) {
            
                auto WSpec = -2.0 * (Wr_fk[fp, q](i, j, i, j) - v_k[q](i, j, i, j)).imag();
                auto num   = bose2(fp * beta) + fermi2(ekq(l) * beta);
                auto den   = f + idelta + fp - ekq(l);
                
//                std::cout << "num: " << num << "\n";
//                std::cout << "den: " << 1.0/den << "\n";
//                std::cout << "mat: " << Ukq(l, i) * dagger(Ukq)(j, l) << "\n";
      
                sigma_fk_dyn[f,k](j,i) += Ukq(l, i) * dagger(Ukq)(j, l) * WSpec * num / den;
                
              }
          
            }
        
          }
          
        }
      }
      
    }
  } 
  
  sigma_fk_dyn *= fmesh.delta() / (kmesh.size() * 2.0 * 3.141592653589793);
  sigma_fk_sta *= 1.0 / kmesh.size();
  
  for (auto const &f : fmesh) {
    for (auto const &k : kmesh) {
      for (int i : range(nb)) {
        for (int j : range(nb)) {
            sigma_fk[f,k](i,j) = sigma_fk_dyn[f,k](i,j);
            //std::cout << sigma_fk_dyn[f,k](i,j) << " - " << sigma_fk_sta[f,k](i,j) << "\n";
        }
      }
    }
  }

  return sigma_fk;
  
}

// ----------------------------------------------------
// gw_sigma_tr
  
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
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[t, r] = arr(idx);

    //for (const auto &[t, r] : g_tr.mesh()) {

    for (const auto &[a, b, c, d] : Wr_tr.target_indices()) {
      sigma_tr[t, r](a, b) += Wr_tr[t, r](a, b, c, d) * g_tr[t, r](c, d);
    }
  }

  return sigma_tr;
}

// ----------------------------------------------------
// gw_sigma_wk_serial_fft

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

} // namespace triqs_tprf
