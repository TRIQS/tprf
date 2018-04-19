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

gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh) {

  auto I = make_unit_matrix<ek_vt::scalar_t>(ek.target_shape()[0]);
  gk_iw_t g0k = make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());
  
  //for (auto const &k : ek.mesh()) {

#pragma omp parallel for 
  for (int idx = 0; idx < ek.mesh().size(); idx++) {
    auto iter = ek.mesh().begin(); iter += idx; auto k = *iter;
    
    for (auto const &w : mesh) g0k[w, k] = inverse((w + mu)*I - ek(k));
  }

  return g0k;
}

gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma) {

  auto mesh = sigma.mesh();
  auto I = make_unit_matrix<ek_vt::scalar_t>(ek.target_shape()[0]);
  gk_iw_t gk =
      make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

  //for (auto const &k : ek.mesh()) {

#pragma omp parallel for 
  for (int idx = 0; idx < ek.mesh().size(); idx++) {
    auto iter = ek.mesh().begin(); iter += idx; auto k = *iter;
    
    for (auto const &w : mesh) gk[w, k] = inverse((w + mu)*I - ek(k) - sigma[w]);
  }

  return gk;
}

  /*
gr_iw_t gr_from_gk(gk_iw_vt gk) {

  auto wmesh = std::get<0>(gk.mesh());
  auto kmesh = std::get<1>(gk.mesh());
  auto lmesh = gf_mesh<cyclic_lattice>{kmesh.domain().lattice(), kmesh.periodization_matrix};

  gr_iw_t gr = make_gf<gr_iw_t::mesh_t::var_t>({wmesh, lmesh}, gk.target());

  for (auto const &w : wmesh) {
    auto _ = var_t{};
    gr[w, _] = inverse_fourier(gk[w, _]);
  }

  return gr;
}

gk_iw_t gk_from_gr(gr_iw_vt gr) {

  auto wmesh = std::get<0>(gr.mesh());
  auto lmesh = std::get<1>(gr.mesh());
  auto kmesh = gf_mesh<brillouin_zone>{brillouin_zone{lmesh.domain()}, lmesh.periodization_matrix};
  
  gk_iw_t gk = make_gf<gk_iw_t::mesh_t::var_t>({wmesh, kmesh}, gr.target());

  for (auto const &w : wmesh) {
    auto _ = var_t{};
    gk[w, _] = fourier(gr[w, _]);
  }

  return gk;
}
 */

gr_iw_t gr_from_gk(gk_iw_vt gwk) {

  auto _ = var_t{};
  auto target = gwk.target();

  auto wmesh = std::get<0>(gwk.mesh());
  auto kmesh = std::get<1>(gwk.mesh());
  auto lmesh = gf_mesh<cyclic_lattice>{kmesh.domain().lattice(), kmesh.periodization_matrix};

  gr_iw_t gwr = make_gf<gr_iw_t::mesh_t::var_t>({wmesh, lmesh}, target);

  auto w0 = *wmesh.begin();
  void * p = _fourier_impl_plan(gwr[w0, _], gwk[w0, _]);
  
#pragma omp parallel for 
  for (int idx = 0; idx < wmesh.size(); idx++) {
    auto iter = wmesh.begin(); iter += idx; auto w = *iter;

    auto gr = make_gf<cyclic_lattice>(lmesh, target);
    auto gk = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
    gk = gwk[w, _];

    _fourier_impl(gr, gk, p);

#pragma omp critical
    gwr[w, _] = gr;

  }

  _fourier_impl_destroy_plan(p);

  return gwr;
}

gk_iw_t gk_from_gr(gr_iw_vt gwr) {

  auto _ = var_t{};
  auto target = gwr.target();

  auto wmesh = std::get<0>(gwr.mesh());
  auto lmesh = std::get<1>(gwr.mesh());
  auto kmesh = gf_mesh<brillouin_zone>{brillouin_zone{lmesh.domain()}, lmesh.periodization_matrix};
  
  gk_iw_t gwk = make_gf<gk_iw_t::mesh_t::var_t>({wmesh, kmesh}, target);

  auto w0 = *wmesh.begin();
  void * p = _fourier_impl_plan(gwk[w0, _], gwr[w0, _]);

#pragma omp parallel for 
  for (int idx = 0; idx < wmesh.size(); idx++) {
    auto iter = wmesh.begin(); iter += idx; auto w = *iter;

    auto gr = make_gf<cyclic_lattice>(lmesh, target);
    auto gk = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
    gr = gwr[w, _];

    _fourier_impl(gk, gr, p);

#pragma omp critical
    gwk[w, _] = gk;

  }

  _fourier_impl_destroy_plan(p);
  
  return gwk;
}

// ----------------------------------------------------
// Transformations: Matsubara frequency <-> imaginary time

/** Set the two lowest order tail coefficients

    for a single-particle Green's function.
    using the 1/w first order coefficient
    and fitting the second order coefficient
    from the value at the lowest Matsubara frequency.
 */

g_iw_vt set_gf_tail(g_iw_vt gw, double s_1 = 1.) {

  auto s = gw.singularity();

  s.zero();

  s(1) = s_1; // set order 1 to the unit matrix

  // get approx of 2nd order from lowest G(w) value

  for (auto const &w : gw.mesh()) {
    s(2) = w * w * (gw[w] - s(1) / w);
    break;
  }

  return gw;
}

gr_tau_t grt_from_grw(gr_iw_vt grw, int ntau) {

  auto wmesh = std::get<0>(grw.mesh());
  auto rmesh = std::get<1>(grw.mesh());

  double beta = wmesh.domain().beta;

  int nw = wmesh.last_index() + 1;
  if( ntau <= 0 ) ntau = 4 * nw;

  gr_tau_t grt = make_gf<gr_tau_t::mesh_t::var_t>(
      {{beta, Fermion, ntau}, rmesh}, grw.target());

  auto _ = var_t{};

  void * p;
  {
    auto r0 = *rmesh.begin();
    auto gw_ab = slice_target_to_scalar(grw[_, r0], 0, 0);
    auto gt_ab = slice_target_to_scalar(grt[_, r0], 0, 0);
    p = _fourier_impl_plan(gt_ab, gw_ab);
  }
  
  //  for (auto const &r : rmesh) {

#pragma omp parallel for 
  for (int idx = 0; idx < rmesh.size(); idx++) {
    auto iter = rmesh.begin(); iter += idx; auto r = *iter;

    double s_1 = 0.;
    if (r.linear_index() == 0)
      s_1 = 1.; // only r=0 has a fermi jump and a 1/i\omega_n tail

    auto gw = make_gf<imfreq>({beta, Fermion, nw}, grw.target());
    auto gt = make_gf<imtime>({beta, Fermion, ntau}, grw.target());

#pragma omp critical
    gw = set_gf_tail(grw[_, r], s_1);

    // ----------------------------------------------------
    for(int a = 0; a < gw.target_shape()[0]; a++) {
      for(int b = 0; b < gw.target_shape()[1]; b++) {
	auto gw_ab = slice_target_to_scalar(gw, a, b);
	auto gt_ab = slice_target_to_scalar(gt, a, b);
	_fourier_impl(gt_ab, gw_ab, p);
      }
    }
    // ----------------------------------------------------
     
    // ----------------------------------------------------
    // serial replacement
    
    //grt[_, r] = inverse_fourier(
    //    set_gf_tail(grw[_, r], s_1)); // Using on the fly tail coeff FIXME
    // ----------------------------------------------------

    // grt[_, r] = inverse_fourier(set_gf_tail(grw[_, r])); // Using on the fly
    // tail coeff FIXME

#pragma omp critical
    grt[_, r] = gt;
  }

  _fourier_impl_destroy_plan(p);

  return grt;
}

} // namespace tprf
