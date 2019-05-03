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

#include "../fourier/fourier.hpp"

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

  g0_wk = mpi::all_reduce(g0_wk);
  return g0_wk;
}

#else
  
g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, gf_mesh<imfreq> mesh) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g0_wk({mesh, e_k.mesh()}, e_k.target_shape());
  
  for (auto const &[w, k] : mpi_view(g0_wk.mesh()))
      g0_wk[w, k] = inverse((w + mu)*I - e_k(k));

  g0_wk = mpi::all_reduce(g0_wk);
  return g0_wk;
}
  
#endif

// ----------------------------------------------------

g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_wk_cvt sigma_wk) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  auto g_wk = make_gf(sigma_wk);

  for (auto const &[w, k] : mpi_view(g_wk.mesh()) ) 
    g_wk[w, k] = inverse((w + mu)*I - e_k[k] - sigma_wk[w, k]);

  g_wk = mpi::all_reduce(g_wk);
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

  g_wk = mpi::all_reduce(g_wk);
  return g_wk;
}

#else
  
g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w) {

  auto mesh = sigma_w.mesh();
  auto I = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g_wk({mesh, e_k.mesh()}, e_k.target_shape());

  for (auto const &[w, k] : mpi_view(g_wk.mesh()) ) 
    g_wk[w, k] = inverse((w + mu)*I - e_k(k) - sigma_w[w]);

  g_wk = mpi::all_reduce(g_wk);
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

  g_w = mpi::all_reduce(g_w);
  g_w /= kmesh.size();

  return g_w;
}

// ----------------------------------------------------
  
#ifdef TPRF_OMP

g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk) {

  auto _ = all_t{};
  auto target = g_wk.target();

  //const auto & [ wmesh, kmesh ] = g_wk.mesh();
  auto wmesh = std::get<0>(g_wk.mesh());
  auto kmesh = std::get<1>(g_wk.mesh());
  auto rmesh = make_adjoint_mesh(kmesh);

  g_wr_t g_wr({wmesh, rmesh}, g_wk.target_shape());

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wk[w0, _]), gf_view(g_wr[w0, _]));

  auto arr = mpi_view(wmesh);

  #pragma omp parallel for 
  for (int idx = 0; idx < arr.size(); idx++) {
    auto & w = arr(idx);

    auto g_r = make_gf<cyclic_lattice>(rmesh, target);
    auto g_k = make_gf<brillouin_zone>(kmesh, target);

    #pragma omp critical
    g_k = g_wk[w, _];

    _fourier_with_plan<0>(gf_const_view(g_k), gf_view(g_r), p);

    #pragma omp critical
    g_wr[w, _] = g_r;

  }

  g_wr = mpi::all_reduce(g_wr);
  return g_wr;
}

#else

g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk) {

  auto [wmesh, kmesh] = g_wk.mesh();
  auto rmesh = make_adjoint_mesh(kmesh);

  g_wr_t g_wr({wmesh, rmesh}, g_wk.target_shape());

  auto _ = all_t{};
  for ( auto const &w : mpi_view(wmesh) )
    g_wr[w, _]() = triqs::gfs::fourier(g_wk[w, _]);  

  g_wr = mpi::all_reduce(g_wr);
  
  return g_wr;
}
  
#endif

// ----------------------------------------------------

#ifdef TPRF_OMP
  
g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr) {

  auto _ = all_t{};
  auto target = g_wr.target();

  //auto [wmesh, rmesh] = g_wr.mesh();
  auto wmesh = std::get<0>(g_wr.mesh());
  auto rmesh = std::get<1>(g_wr.mesh());
  auto kmesh = make_adjoint_mesh(rmesh);
  
  g_wk_t g_wk({wmesh, kmesh}, g_wr.target_shape());

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wr[w0, _]), gf_view(g_wk[w0, _]));

  auto arr = mpi_view(wmesh);

  #pragma omp parallel for 
  for (int idx = 0; idx < arr.size(); idx++) {
    auto & w = arr(idx);

    auto g_r = make_gf<cyclic_lattice>(rmesh, target);
    auto g_k = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
    g_r = g_wr[w, _];

    _fourier_with_plan<0>(gf_const_view(g_r), gf_view(g_k), p);

#pragma omp critical
    g_wk[w, _] = g_k;

  }
  
  g_wk = mpi::all_reduce(g_wk);
  return g_wk;
}

#else
  
g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr) {

  auto [wmesh, rmesh] = g_wr.mesh();
  auto kmesh = make_adjoint_mesh(rmesh);
  
  g_wk_t g_wk({wmesh, kmesh}, g_wr.target_shape());

  auto _ = all_t{};
  for (auto const &w : mpi_view(wmesh))
    g_wk[w, _]() = triqs::gfs::fourier(g_wr[w, _]);

  g_wk = mpi::all_reduce(g_wk);
  
  return g_wk;
}

#endif
  
// ----------------------------------------------------
// Transformations: Matsubara frequency <-> imaginary time

g_wr_t fourier_tr_to_wr(g_tr_cvt g_tr, int nw) {
  std::cout << "WARNING: fourier_tr_to_wr is not parallellized. FIXME\n";

  auto tmesh = std::get<0>(g_tr.mesh());
  double beta = tmesh.domain().beta;
  auto S = tmesh.domain().statistic;

  if( nw <= 0 ) nw = tmesh.size() / 4;

  auto wmesh = gf_mesh<imfreq>(beta, S, nw);
  
  auto g_wr = make_gf_from_fourier<0>(g_tr, wmesh);
  
  return g_wr;
}

#ifdef TPRF_OMP
  
g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt) {

  auto wmesh = std::get<0>(g_wr.mesh());
  auto rmesh = std::get<1>(g_wr.mesh());

  double beta = wmesh.domain().beta;
  auto S = wmesh.domain().statistic;

  int nw = wmesh.last_index() + 1;
  if( nt <= 0 ) nt = 4 * nw;

  g_tr_t g_tr({{beta, S, nt}, rmesh}, g_wr.target_shape());

  auto tmesh = std::get<0>(g_tr.mesh());

  auto _ = all_t{};

  auto r0 = *rmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(g_wr[_, r0]), gf_view(g_tr[_, r0]));

  auto arr = mpi_view(rmesh);

#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto & r = arr(idx);

    auto g_w = make_gf<imfreq>(wmesh, g_wr.target());
    auto g_t = make_gf<imtime>(tmesh, g_wr.target());

#pragma omp critical
    g_w = g_wr[_, r];

    _fourier_with_plan<0>(gf_const_view(g_w), gf_view(g_t), p);

#pragma omp critical
    g_tr[_, r] = g_t;

  }

  g_tr = mpi::all_reduce(g_tr);
  return g_tr;
}

#else
  
g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt) {

  auto wmesh = std::get<0>(g_wr.mesh());
  auto rmesh = std::get<1>(g_wr.mesh());

  double beta = wmesh.domain().beta;

  int nw = wmesh.last_index() + 1;
  if( nt <= 0 ) nt = 4 * nw;

  g_tr_t g_tr({{beta, wmesh.domain().statistic, nt}, rmesh}, g_wr.target_shape());

  auto _ = all_t{};
  auto r0 = *rmesh.begin();
  auto zero_tail = make_zero_tail(g_wr[_, r0], 2);
  auto zero_tail_r0 = make_zero_tail(g_wr[_, r0], 2);

  zero_tail_r0(1, range(), range()) =
    make_unit_matrix<dcomplex>(g_wr.target_shape()[0]);
  
  for (auto const &r : rmesh) {
    if(r.linear_index() == 0)
      g_tr[_, r]() = triqs::gfs::fourier<0>(g_wr[_, r], make_const_view(zero_tail_r0));
    else
      g_tr[_, r]() = triqs::gfs::fourier<0>(g_wr[_, r], make_const_view(zero_tail));
  }

  return g_tr;
}

#endif
  
} // namespace triqs_tprf
