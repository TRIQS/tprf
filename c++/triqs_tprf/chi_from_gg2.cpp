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

#include "chi_from_gg2.hpp"

using namespace nda::clef;

namespace {
placeholder<0> a;
placeholder<1> b;
placeholder<2> c;
placeholder<3> d;

placeholder<4> Omega;
placeholder<5> n;
placeholder<6> np;
} // namespace

namespace triqs_tprf {

// ----------------------------------------------------
// Particle-hole (PH)

template <> g2_iw_t chi0_from_gg2<Channel_t::PH>(g_iw_cvt g, g2_iw_cvt g2) {
  double beta = g.mesh().domain().beta;
  auto chi0 = make_gf(g2.mesh(), g2.target(), g2.memory_layout());

  chi0(Omega, n, np)(a, b, c, d)
      << -beta * kronecker(n, np) * g(n)(d, a) * g(Omega + n)(b, c);

  return chi0;
}

template <> g2_iw_t chi_from_gg2<Channel_t::PH>(g_iw_cvt g, g2_iw_cvt g2) {
  double beta = g.mesh().domain().beta;
  auto chi = make_gf(g2.mesh(), g2.target(), g2.memory_layout());

  chi(Omega, n, np)(a, b, c, d)
      << g2(Omega, n, np)(a, b, c, d) -
             beta * kronecker(Omega) * g(n)(b, a) * g(np)(d, c);

  return chi;
}

// ----------------------------------------------------
// Particle-particle (PP)

template <> g2_iw_t chi0_from_gg2<Channel_t::PP>(g_iw_cvt g, g2_iw_cvt g2) {
  double beta = g.mesh().domain().beta;
  auto chi0 = make_gf(g2.mesh(), g2.target(), g2.memory_layout());

  chi0(Omega, n, np)(a, b, c, d)
      << -beta * kronecker(n, np) * g(n)(d, a) * g(Omega - n)(b, c);

  return chi0;
}

template <> g2_iw_t chi_from_gg2<Channel_t::PP>(g_iw_cvt g, g2_iw_cvt g2) {
  double beta = g.mesh().domain().beta;
  auto chi = make_gf(g2.mesh(), g2.target(), g2.memory_layout());

  chi(Omega, n, np)(a, b, c, d)
      << g2(Omega, n, np)(a, b, c, d) -
             beta * kronecker(n + np, Omega) * g(n)(b, a) * g(np)(d, c);

  return chi;
}


// ----------------------------------------------------
// functions for (easier) Python wrapping
  
g2_iw_t chi0_from_gg2_PH(g_iw_vt g, g2_iw_vt g2) { return chi0_from_gg2<Channel_t::PH>(g, g2); }
g2_iw_t chi0_from_gg2_PP(g_iw_vt g, g2_iw_vt g2) { return chi0_from_gg2<Channel_t::PP>(g, g2); }

g2_iw_t chi_from_gg2_PH(g_iw_vt g, g2_iw_vt g2) { return chi_from_gg2<Channel_t::PH>(g, g2); }
g2_iw_t chi_from_gg2_PP(g_iw_vt g, g2_iw_vt g2) { return chi_from_gg2<Channel_t::PP>(g, g2); }

} // namespace triqs_tprf
