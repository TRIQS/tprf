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

#include "freq_conv.hpp"

#include <triqs/clef.hpp>
using namespace triqs::clef;

namespace {
placeholder<0> a;
placeholder<1> b;
placeholder<2> c;
placeholder<3> d;

placeholder<4> w;
placeholder<5> n;
placeholder<6> np;

placeholder<7> n1;
placeholder<8> n2;
placeholder<9> n3;
} // namespace

namespace triqs_tprf {
  
// ----------------------------------------------------
// AB matrix_valued from block_gf

g_iw_t block_iw_AB_to_matrix_valued(b_g_iw_vt bg_AB) {

  int n = 0;
  for( auto bidx : range(bg_AB.size()) ) {
    auto g_AB = bg_AB(bidx);
    n += g_AB.target_shape()[0];
  }

  gf<imfreq, matrix_valued> g{bg_AB(0).mesh(), {n, n}};
  g *= 0.0;

  int idx = 0;
  for( auto bidx : range(bg_AB.size()) ) {
    auto g_AB = bg_AB(bidx);
    int size = g_AB.target_shape()[0];
    
    g.data()(range(), range(idx, idx+size), range(idx, idx+size))
      = g_AB.data()(range(), range(), range());

    idx += size;
  }
  
  return g;
}
  
// ----------------------------------------------------
// AABB to ABBA

/**

# ABBA from AABB

General expressions

$$
G(\bar{a}a\bar{b}b) = 
\langle \bar{a} a \bar{b} b \rangle =
- \langle \bar{a} b \bar{b} a \rangle =
- G(\bar{a} b \bar{b} a)
$$

$$
G(\bar{a}b\bar{c}d) = -G(\bar{a}d\bar{c}b)
$$

$$
G_{\bar{a}b\bar{c}d}(\nu_1,\nu_2,\nu_3,\nu_1-\nu_2+\nu_3) = 
-G(\bar{a}d\bar{c}b)(\nu_1, \nu_1-\nu_2+\nu_3, \nu_3, \nu_2)
$$

 */
  
void block_3nu_AABB_to_tensor_valued(b_g2_iw_vt bg2_AABB, g2_iw_vt g2) {
  
  g2 *= 0.0;

  // set AABB components

  // Clef with repeated indices is broken: https://github.com/TRIQS/triqs/issues/475
  //g2(n1, n2, n3)(b1, b1, b2, b2) << bg2_AABB(b1, b2)(n1, n2, n3)(0, 0, 0, 0); // set AABB comp

  // kronecker functionc fix clef deficiency
  g2(n1, n2, n3)(a, b, c, d) << kronecker(a,b) * kronecker(c,d) *
    bg2_AABB(a, c)(n1, n2, n3)(0, 0, 0, 0); 

  // ABBA from AABB

  // Clef with repeated indices is broken: https://github.com/TRIQS/triqs/issues/475
  //g2(n1, n2, n3)(b1, b2, b2, b1) << (1-kronecker(b1, b2)) * -g2(n1, n1 - n2 + n3, n3)(b1, b1, b2, b2);

  g2(n1, n2, n3)(a, b, c, d) << g2(n1, n2, n3)(a, b, c, d) +
    kronecker(a,d) * kronecker(b,c) * (1 - kronecker(a, c)) *
    -g2(n1, n1 - n2 + n3, n3)(a, d, c, b);
  
}

// ----------------------------------------------------

void get_magnetic_component(g2_iw_vt g2, g2_iw_vt g2_m) {
  g2_m(n1, n2, n3) << g2(n1, n2, n3)(0, 0, 0, 0) - g2(n1, n2, n3)(0, 0, 1, 1);
}
  
// ----------------------------------------------------
// Particle-hole (PH)

template <> void from_3nu<Channel_t::PH>(g2_iw_vt g2_ch, g2_iw_cvt g2) {
  g2_ch(w, n, np) << g2(n, n + w, np + w);
}

// ----------------------------------------------------
// Particle-hole-bar (PHbar)

template <> void from_3nu<Channel_t::PH_bar>(g2_iw_vt g2_ch, g2_iw_cvt g2) {
  g2_ch(w, n, np) << g2(n, np, np + w);
}

// ----------------------------------------------------
// Particle-particle (PP)

template <> void from_3nu<Channel_t::PP>(g2_iw_vt g2_ch, g2_iw_cvt g2) {
  g2_ch(w, n, np) << g2(n, np, w - np);
}

void from_3nu_PH(g2_iw_vt g2_ch, g2_iw_vt g2) { from_3nu<Channel_t::PH>(g2_ch, g2); }
void from_3nu_PH_bar(g2_iw_vt g2_ch, g2_iw_vt g2) { from_3nu<Channel_t::PH_bar>(g2_ch, g2); }
void from_3nu_PP(g2_iw_vt g2_ch, g2_iw_vt g2) { from_3nu<Channel_t::PP>(g2_ch, g2); }

} // namespace triqs_tprf
