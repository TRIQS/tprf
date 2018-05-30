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

#include "common.hpp"
#include "gf.hpp"

namespace tprf {
  
// ----------------------------------------------------
// g

#ifdef TPRF_OMP

gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(ek.target_shape()[0]);
  gk_iw_t g0k = make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());
  
#pragma omp parallel for 
  for (int idx = 0; idx < ek.mesh().size(); idx++) {
    auto iter = ek.mesh().begin(); iter += idx; auto k = *iter;
    
    for (auto const &w : mesh) g0k[w, k] = inverse((w + mu)*I - ek(k));
  }

  return g0k;
}

#else
  
gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(ek.target_shape()[0]);
  gk_iw_t g0k = make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());
  
  for (auto const &k : ek.mesh())
    for (auto const &w : mesh) g0k[w, k] = inverse((w + mu)*I - ek(k));

  return g0k;
}

#endif

// ----------------------------------------------------
  
#ifdef TPRF_OMP

gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma) {

  auto mesh = sigma.mesh();
  auto I = make_unit_matrix<ek_vt::scalar_t>(ek.target_shape()[0]);
  gk_iw_t gk =
      make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

#pragma omp parallel for 
  for (int idx = 0; idx < ek.mesh().size(); idx++) {
    auto iter = ek.mesh().begin(); iter += idx; auto k = *iter;
    
    for (auto const &w : mesh) gk[w, k] = inverse((w + mu)*I - ek(k) - sigma[w]);
  }

  return gk;
}

#else
  
gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma) {

  auto mesh = sigma.mesh();
  auto I = make_unit_matrix<ek_vt::scalar_t>(ek.target_shape()[0]);
  gk_iw_t gk =
      make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

  for (auto const &k : ek.mesh())
    for (auto const &w : mesh) gk[w, k] = inverse((w + mu)*I - ek(k) - sigma[w]);

  return gk;
}

#endif

// ----------------------------------------------------
  
#ifdef TPRF_OMP

gr_iw_t gr_from_gk(gk_iw_vt gwk) {

  auto _ = all_t{};
  auto target = gwk.target();

  //const auto & [ wmesh, kmesh ] = gwk.mesh();
  auto wmesh = std::get<0>(gwk.mesh());
  auto kmesh = std::get<1>(gwk.mesh());
  auto rmesh = make_adjoint_mesh(kmesh);

  gr_iw_t gwr = make_gf<gr_iw_t::mesh_t::var_t>({wmesh, rmesh}, target);

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(gwk[w0, _]), gf_view(gwr[w0, _]));

#pragma omp parallel for 
  for (int idx = 0; idx < wmesh.size(); idx++) {
    auto iter = wmesh.begin(); iter += idx; auto w = *iter;

    auto gr = make_gf<cyclic_lattice>(rmesh, target);
    auto gk = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
    gk = gwk[w, _];

    _fourier_with_plan<0>(gf_const_view(gk), gf_view(gr), p);

#pragma omp critical
    gwr[w, _] = gr;

  }

  return gwr;
}

#else
  
gr_iw_t gr_from_gk(gk_iw_vt gwk) {

  auto [wmesh, kmesh] = gwk.mesh();
  auto rmesh = make_adjoint_mesh(kmesh);

  gr_iw_t gwr = make_gf<gr_iw_t::mesh_t::var_t>({wmesh, rmesh}, gwk.target());

  auto _ = all_t{};
  for (auto const &w : wmesh) gwr[w, _]() = fourier(gwk[w, _]);

  return gwr;
}

#endif

// ----------------------------------------------------

#ifdef TPRF_OMP
  
gk_iw_t gk_from_gr(gr_iw_vt gwr) {

  auto _ = all_t{};
  auto target = gwr.target();

  //auto [wmesh, rmesh] = gwr.mesh();
  auto wmesh = std::get<0>(gwr.mesh());
  auto rmesh = std::get<1>(gwr.mesh());
  auto kmesh = make_adjoint_mesh(rmesh);
  
  gk_iw_t gwk = make_gf<gk_iw_t::mesh_t::var_t>({wmesh, kmesh}, target);

  auto w0 = *wmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(gwr[w0, _]), gf_view(gwk[w0, _]));

#pragma omp parallel for 
  for (int idx = 0; idx < wmesh.size(); idx++) {
    auto iter = wmesh.begin(); iter += idx; auto w = *iter;

    auto gr = make_gf<cyclic_lattice>(rmesh, target);
    auto gk = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
    gr = gwr[w, _];

    _fourier_with_plan<0>(gf_const_view(gr), gf_view(gk), p);

#pragma omp critical
    gwk[w, _] = gk;

  }
  
  return gwk;
}

#else
  
gk_iw_t gk_from_gr(gr_iw_vt gwr) {

  auto [wmesh, rmesh] = gwr.mesh();
  auto kmesh = make_adjoint_mesh(rmesh);
  
  gk_iw_t gwk = make_gf<gk_iw_t::mesh_t::var_t>({wmesh, kmesh}, gwr.target());

  auto _ = all_t{};
  for (auto const &w : wmesh) gwk[w, _]() = fourier(gwr[w, _]);
  
  return gwk;
}

#endif
  
// ----------------------------------------------------
// Transformations: Matsubara frequency <-> imaginary time

#ifdef TPRF_OMP
  
gr_tau_t grt_from_grw(gr_iw_vt grw, int ntau) {

  auto wmesh = std::get<0>(grw.mesh());
  auto rmesh = std::get<1>(grw.mesh());

  double beta = wmesh.domain().beta;
  auto S = wmesh.domain().statistic;

  int nw = wmesh.last_index() + 1;
  if( ntau <= 0 ) ntau = 4 * nw;

  gr_tau_t grt = make_gf<gr_tau_t::mesh_t::var_t>(
      {{beta, S, ntau}, rmesh}, grw.target());

  auto tmesh = std::get<0>(grt.mesh());

  auto _ = all_t{};

  auto r0 = *rmesh.begin();
  auto p = _fourier_plan<0>(gf_const_view(grw[_, r0]), gf_view(grt[_, r0]));

#pragma omp parallel for 
  for (int idx = 0; idx < rmesh.size(); idx++) {
    auto iter = rmesh.begin(); iter += idx; auto r = *iter;

    auto gw = make_gf<imfreq>(wmesh, grw.target());
    auto gt = make_gf<imtime>(tmesh, grw.target());

#pragma omp critical
    gw = grw[_, r];

    _fourier_with_plan<0>(gf_const_view(gw), gf_view(gt), p);

#pragma omp critical
    grt[_, r] = gt;

  }

  return grt;
}

#else
  
gr_tau_t grt_from_grw(gr_iw_vt grw, int ntau) {

  auto wmesh = std::get<0>(grw.mesh());
  auto rmesh = std::get<1>(grw.mesh());

  double beta = wmesh.domain().beta;

  int nw = wmesh.last_index() + 1;
  if( ntau <= 0 ) ntau = 4 * nw;

  gr_tau_t grt = make_gf<gr_tau_t::mesh_t::var_t>(
    {{beta, wmesh.domain().statistic, ntau}, rmesh}, grw.target());

  auto _ = all_t{};
  for (auto const &r : rmesh) grt[_, r]() = fourier<0>(grw[_, r]);

  return grt;
}

#endif
  
} // namespace tprf
