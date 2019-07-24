 /*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
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

#include <triqs/arrays/linalg/det_and_inverse.hpp>
using triqs::arrays::inverse;

#include "../mpi.hpp"
#include "common.hpp"
#include "gf.hpp"

#include "fourier_gf.hpp"

namespace triqs_tprf {

  namespace {
    using namespace fourier;
  }
  
// ----------------------------------------------------
// g

#ifdef TPRF_OMP

g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, gf_mesh<imfreq> mesh) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g0_wk({mesh, e_k.mesh()}, e_k.target_shape());
    
  auto arr = mpi_view(g0_wk.mesh());

#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);
    g0_wk[w, k] = inverse((w + mu)*I - e_k(k));      
  }

  g0_wk = mpi_all_reduce(g0_wk);  
  return g0_wk;
}

#else
  
g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, gf_mesh<imfreq> mesh) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g0_wk({mesh, e_k.mesh()}, e_k.target_shape());
  
  for (auto const &[w, k] : mpi_view(g0_wk.mesh()))
      g0_wk[w, k] = inverse((w + mu)*I - e_k(k));

  g0_wk = mpi_all_reduce(g0_wk);
  return g0_wk;
}
  
#endif

// ----------------------------------------------------

g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_wk_cvt sigma_wk) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  auto g_wk = make_gf(sigma_wk);

  for (auto const &[w, k] : mpi_view(g_wk.mesh()) ) 
    g_wk[w, k] = inverse((w + mu)*I - e_k[k] - sigma_wk[w, k]);

  g_wk = mpi_all_reduce(g_wk);
  return g_wk;
}
  
#ifdef TPRF_OMP

g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w) {

  auto mesh = sigma_w.mesh();
  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g_wk({sigma_w.mesh(), e_k.mesh()}, e_k.target_shape());

  auto arr = mpi_view(g_wk.mesh());

#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);
    g_wk[w, k] = inverse((w + mu)*I - e_k(k) - sigma_w[w]);
  }

  g_wk = mpi_all_reduce(g_wk);  
  return g_wk;
}

#else
  
g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w) {

  auto mesh = sigma_w.mesh();
  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g_wk({mesh, e_k.mesh()}, e_k.target_shape());

  for (auto const &[w, k] : mpi_view(g_wk.mesh()) ) 
    g_wk[w, k] = inverse((w + mu)*I - e_k(k) - sigma_w[w]);

  g_wk = mpi_all_reduce(g_wk);
  return g_wk;
}

#endif

g_w_t lattice_dyson_g_w(double mu, e_k_cvt e_k, g_w_cvt sigma_w) {

  auto wmesh = sigma_w.mesh();
  auto kmesh = e_k.mesh();

  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_w_t g_w(wmesh, e_k.target_shape());

  auto wkmesh = gk_iw_t::mesh_t{{wmesh, kmesh}};
  
  for (auto const &[w, k] : mpi_view(wkmesh) ) 
    g_w[w] += inverse((w + mu)*I - e_k(k) - sigma_w[w]);

  g_w = mpi_all_reduce(g_w);
  g_w /= kmesh.size();

  return g_w;
}

// ----------------------------------------------------
// Transformations: real space <-> reciprocal space 
  
g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk) {
  auto g_wr = fourier_wk_to_wr_general_target(g_wk);
  return g_wr;
}

g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr) {
  auto g_wk = fourier_wr_to_wk_general_target(g_wr);
  return g_wk;
}

// ----------------------------------------------------
// Transformations: Matsubara frequency <-> imaginary time

g_wr_t fourier_tr_to_wr(g_tr_cvt g_tr, int nw) {
  auto g_wr = fourier_tr_to_wr_general_target(g_tr, nw);
  return g_wr;
}

g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt) {
  auto g_tr = fourier_wr_to_tr_general_target(g_wr, nt);
  return g_tr;
}

} // namespace triqs_tprf
