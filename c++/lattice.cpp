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

#include "lattice.hpp"

#include <triqs/clef.hpp>
using namespace triqs::clef;

namespace {
placeholder<0> iw;
placeholder<1> inu;
placeholder<2> k;
placeholder<3> r;

placeholder<4> a;
placeholder<5> b;
placeholder<6> c;
placeholder<7> d;
} // namespace

#include <triqs/arrays/linalg/det_and_inverse.hpp>
using triqs::arrays::inverse;

namespace tprf {

  void test_function() {
    std::cout << "Hello lattice module!\n";
  }

  void test_function_tb(triqs::lattice::tight_binding tb) {
    std::cout << "Hello lattice module tb arg!\n";
  }

  void test_function_bz(triqs::lattice::brillouin_zone bz) {
    std::cout << "Hello lattice module bz arg!\n";
  }

  void test_function_ek(ek_vt ek){
    std::cout << "Hello lattice module ek arg!\n";
  }

  //void test_function_gk(gk_iw_cvt gk) { std::cout << "Hello lattice module gk arg!\n"; }

  //void test_function_gk(gf_const_view<cartesian_product<imfreq, brillouin_zone>> gk)  { std::cout << "Hello lattice module gk arg!\n"; }
  //void test_function_gk(gf_view<cartesian_product<imfreq, brillouin_zone>> gk)  { std::cout << "Hello lattice module gk arg!\n"; }

  void test_function_gk(gk_iw_vt gk)  { std::cout << "Hello lattice module gk arg!\n"; }

  
  
gk_iw_t g0k_from_ek(double mu, ek_cvt ek, g_iw_t::mesh_t mesh) {
  gk_iw_t g0k = make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

  g0k(inu, k)(a, b) << kronecker(a, b) * (inu + mu) - ek(k)(a, b);

  for (auto const &kidx : std::get<1>(g0k.mesh())) {
    auto _ = var_t{};
    g0k[_][kidx] = inverse(g0k[_][kidx]);
  }

  // g0k = inverse(g0k);  // does not work, see TRIQS issue #463
  return g0k;
}

gk_iw_t gk_from_ek_sigma(double mu, ek_cvt ek, g_iw_cvt sigma) {

  gk_iw_t gk =
      make_gf<gk_iw_t::mesh_t::var_t>({sigma.mesh(), ek.mesh()}, ek.target());

  gk(inu, k)(a, b) << kronecker(a, b) * (inu + mu) - ek(k)(a, b) -
                          sigma(inu)(a, b);

  for (auto const &kidx : std::get<1>(gk.mesh())) {
    auto _ = var_t{};
    gk[_][kidx] = inverse(gk[_][kidx]);
  }

  // gk = inverse(gk);  // does not work, see TRIQS issue #463
  return gk;
}

gr_iw_t gr_from_gk(gk_iw_vt gk) {
  int nk = std::get<1>(gk.mesh()).get_dimensions()[0];

  gr_iw_t gr = make_gf<gr_iw_t::mesh_t::var_t>(
      {std::get<0>(gk.mesh()), {nk, nk}}, gk.target());

  auto _ = var_t{};
  for (auto const &nidx : std::get<0>(gr.mesh()))
    gr[nidx][_] = inverse_fourier(gk[nidx][_]);

  return gr;
}

gk_iw_t gk_from_gr(gr_iw_cvt gr, brillouin_zone const &bz) {
  int nk = std::get<1>(gr.mesh()).get_dimensions()[0];

  gk_iw_t gk = make_gf<gk_iw_t::mesh_t::var_t>(
      {std::get<0>(gr.mesh()), {bz, nk}}, gr.target());

  auto _ = var_t{};
  for (auto const &nidx : std::get<0>(gk.mesh()))
    gk[nidx][_] = fourier(gr[nidx][_]);

  return gk;
}

chi0r_t chi0r_from_gr_PH(int nw, int nnu, gr_iw_cvt gr) {

  int nb = gr.target().shape()[0];
  int nk = std::get<1>(gr.mesh()).get_dimensions()[0];
  double beta = std::get<0>(gr.mesh()).domain().beta;

  chi0r_t chi0r{{{beta, Boson, nw}, {beta, Fermion, nnu}, {nk, nk}},
                {nb, nb, nb, nb}};

  chi0r(iw, inu, r)(a, b, c, d) << gr(inu, r)(d, a) * gr(inu + iw, -r)(b, c);

  return chi0r;
}

chi0r_t chi0r_from_chi0q(chi0q_cvt chi0q) {
  auto mb = std::get<0>(chi0q.mesh());
  auto mf = std::get<1>(chi0q.mesh());
  int nk = std::get<2>(chi0q.mesh()).get_dimensions()[0];

  auto chi0r =
      make_gf<chi0r_t::mesh_t::var_t>({mb, mf, {nk, nk}}, chi0q.target());

  for (auto const &widx : std::get<0>(chi0r.mesh())) {
    for (auto const &nidx : std::get<1>(chi0r.mesh())) {
      auto _ = var_t{};
      chi0r[widx][nidx][_] = inverse_fourier(chi0q[widx][nidx][_]);
    }
  }
  return chi0r;
}

chi0q_t chi0q_from_chi0r(chi0r_cvt chi0r, brillouin_zone const &bz) {
  auto mb = std::get<0>(chi0r.mesh());
  auto mf = std::get<1>(chi0r.mesh());
  int nk = std::get<2>(chi0r.mesh()).get_dimensions()[0];

  auto chi0q =
      make_gf<chi0q_t::mesh_t::var_t>({mb, mf, {bz, nk}}, chi0r.target());

  for (auto const &widx : std::get<0>(chi0q.mesh())) {
    for (auto const &nidx : std::get<1>(chi0q.mesh())) {
      auto _ = var_t{};
      chi0q[widx][nidx][_] = fourier(chi0r[widx][nidx][_]);
    }
  }
  return chi0q;
}

} // namespace tprf
