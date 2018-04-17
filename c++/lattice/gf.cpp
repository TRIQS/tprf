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

  gk_iw_t g0k = make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

  matrix<ek_vt::scalar_t> I(ek.target_shape());
  I(a, b) << kronecker(a, b);
  
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
  
  gk_iw_t gk =
      make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

  matrix<ek_vt::scalar_t> I(ek.target_shape());
  I(a, b) << kronecker(a, b);

  //for (auto const &k : ek.mesh()) {

#pragma omp parallel for 
  for (int idx = 0; idx < ek.mesh().size(); idx++) {
    auto iter = ek.mesh().begin(); iter += idx; auto k = *iter;
    
    for (auto const &w : mesh) gk[w, k] = inverse((w + mu)*I - ek(k) - sigma[w]);
  }

  return gk;
}
  
gr_iw_t gr_from_gk(gk_iw_vt gk) {

  auto wmesh = std::get<0>(gk.mesh());
  auto kmesh = std::get<1>(gk.mesh());
  auto lmesh = gf_mesh<cyclic_lattice>{kmesh.domain().lattice(), kmesh.periodization_matrix};

  gr_iw_t gr = make_gf<gr_iw_t::mesh_t::var_t>({wmesh, lmesh}, gk.target());

  for (auto const &w : wmesh) {

  /*
#pragma omp parallel for 
  for (int idx = 0; idx < wmesh.size(); idx++) {
    auto iter = wmesh.begin(); iter += idx; auto w = *iter;
  */

    auto _ = var_t{};
#pragma omp single
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
#pragma omp single
    gk[w, _] = fourier(gr[w, _]);
  }

  return gk;
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

gr_tau_t grt_from_grw(gr_iw_vt grw) {

  auto wmesh = std::get<0>(grw.mesh());
  auto rmesh = std::get<1>(grw.mesh());

  double beta = wmesh.domain().beta;

  int nw = wmesh.last_index() + 1;
  int ntau = 4 * nw;

  gr_tau_t grt = make_gf<gr_tau_t::mesh_t::var_t>(
      {{beta, Fermion, ntau}, rmesh}, grw.target());

  auto _ = var_t{};
  for (auto const &r : rmesh) {

    auto gw = grw[_, r];

    // std::cout << "r = " << r << "\n";
    // std::cout << "ridx = " << r.linear_index() << "\n";

    double s_1 = 0.;
    if (r.linear_index() == 0)
      s_1 = 1.; // only r=0 has a fermi jump and a 1/i\omega_n tail
    // std::cout << "s_1 = " << s_1 << "\n";

    grt[_, r] = inverse_fourier(
        set_gf_tail(grw[_, r], s_1)); // Using on the fly tail coeff FIXME

    // grt[_, r] = inverse_fourier(set_gf_tail(grw[_, r])); // Using on the fly
    // tail coeff FIXME
  }

  return grt;
}

} // namespace tprf
