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

#include "linalg.hpp"

namespace triqs_tprf {

 // ----------------------------------------------------

 template <Channel_t CH> g2_nn_t inverse(g2_nn_cvt g) {

  channel_grouping<CH> chg;
  auto g_inv = make_gf(g.mesh(), g.target(), g.memory_layout());

  auto g_w = make_gf(g, chg.memory_layout());
  auto g_w_inv = make_gf(g_w.mesh(), g_w.target(), g_w.memory_layout());

  auto mat = chg.matrix_view(g_w.data());
  auto mat_inv = chg.matrix_view(g_w_inv.data());

  mat_inv = inverse(mat);

  g_inv = g_w_inv;

  return g_inv;
 }
  
 /// Inverse: [G]^{-1}, Two-particle response-function inversion
 template <Channel_t CH> g2_iw_t inverse(g2_iw_cvt g) {

   //channel_grouping<CH> chg;
  auto g_inv = make_gf(g.mesh(), g.target(), g.memory_layout());

  for (auto const &w : std::get<0>(g.mesh())) {
   auto _ = all_t{};
   g_inv[w, _, _] = inverse<CH>(g[w, _, _]);
  }
  return g_inv;
 }

 // ----------------------------------------------------
  
 template <Channel_t CH> g2_nn_t product(g2_nn_cvt A, g2_nn_cvt B) {

   channel_grouping<CH> chg;
   
   auto C = make_gf(A.mesh(), A.target(), A.memory_layout());

   auto A_w = make_gf(A, chg.memory_layout());
   auto B_w = make_gf(B, chg.memory_layout());
   auto C_w = make_gf(A_w.mesh(), A_w.target(), A_w.memory_layout());

   auto A_mat = chg.matrix_view(A_w.data());
   auto B_mat = chg.matrix_view(B_w.data());
   auto C_mat = chg.matrix_view(C_w.data());

   C_mat = A_mat * B_mat;

   C = C_w;

   return C;
 }

 /// product: C = A * B, two-particle response-function product
 template <Channel_t CH> g2_iw_t product(g2_iw_cvt A, g2_iw_cvt B) {

  //channel_grouping<CH> chg;
  // check that A and B are compatible!

  auto C = make_gf(A.mesh(), A.target(), A.memory_layout());

  for (auto const &w : std::get<0>(A.mesh())) {
   auto _ = all_t{};
   C[w, _, _] = product<CH>(A[w, _, _], B[w, _, _]);
  }
  return C;
 }

 // ----------------------------------------------------

 template <Channel_t CH> g2_nn_t identity(g2_nn_cvt g) {

  channel_grouping<CH> chg;
  auto I = make_gf(g.mesh(), g.target(), g.memory_layout());
  auto I_w = make_gf(I, chg.memory_layout());
  auto I_mat = chg.matrix_view(I_w.data());

  I_mat = 1.0; // This sets a triqs::arrays::matrix to the identity matrix...
  I = I_w;

  return I;
 }

 /// Identity: 1, identity two-particle response-function
 template <Channel_t CH> g2_iw_t identity(g2_iw_cvt g) {

  auto I = make_gf(g.mesh(), g.target(), g.memory_layout());

  for (auto const &w : std::get<0>(g.mesh())) {
   auto _ = all_t{};
   I[w, _, _] = identity<CH>(I[w, _, _]);
  }
  return I;
 }
  
 // ----------------------------------------------------

 template g2_nn_t inverse<Channel_t::PH>(g2_nn_cvt);
 template g2_nn_t inverse<Channel_t::PH_bar>(g2_nn_cvt);
 template g2_nn_t inverse<Channel_t::PP>(g2_nn_cvt);

 g2_nn_t inverse_PH(g2_nn_vt g) { return inverse<Channel_t::PH>(g); }
 g2_nn_t inverse_PP(g2_nn_vt g) { return inverse<Channel_t::PP>(g); }
 g2_nn_t inverse_PH_bar(g2_nn_vt g) { return inverse<Channel_t::PH_bar>(g); }

 template g2_iw_t inverse<Channel_t::PH>(g2_iw_cvt);
 template g2_iw_t inverse<Channel_t::PH_bar>(g2_iw_cvt);
 template g2_iw_t inverse<Channel_t::PP>(g2_iw_cvt);

 g2_iw_t inverse_PH(g2_iw_vt g) { return inverse<Channel_t::PH>(g); }
 g2_iw_t inverse_PP(g2_iw_vt g) { return inverse<Channel_t::PP>(g); }
 g2_iw_t inverse_PH_bar(g2_iw_vt g) { return inverse<Channel_t::PH_bar>(g); }

 // ----------------------------------------------------
  
 template g2_nn_t product<Channel_t::PH>(g2_nn_cvt, g2_nn_cvt);
 template g2_nn_t product<Channel_t::PH_bar>(g2_nn_cvt, g2_nn_cvt);
 template g2_nn_t product<Channel_t::PP>(g2_nn_cvt, g2_nn_cvt);

 g2_nn_t product_PH(g2_nn_vt A, g2_nn_vt B) { return product<Channel_t::PH>(A, B); }
 g2_nn_t product_PP(g2_nn_vt A, g2_nn_vt B) { return product<Channel_t::PP>(A, B); }
 g2_nn_t product_PH_bar(g2_nn_vt A, g2_nn_vt B) { return product<Channel_t::PH_bar>(A, B); }

 template g2_iw_t product<Channel_t::PH>(g2_iw_cvt, g2_iw_cvt);
 template g2_iw_t product<Channel_t::PH_bar>(g2_iw_cvt, g2_iw_cvt);
 template g2_iw_t product<Channel_t::PP>(g2_iw_cvt, g2_iw_cvt);

 g2_iw_t product_PH(g2_iw_vt A, g2_iw_vt B) { return product<Channel_t::PH>(A, B); }
 g2_iw_t product_PP(g2_iw_vt A, g2_iw_vt B) { return product<Channel_t::PP>(A, B); }
 g2_iw_t product_PH_bar(g2_iw_vt A, g2_iw_vt B) { return product<Channel_t::PH_bar>(A, B); }
 
 // ----------------------------------------------------

 template g2_nn_t identity<Channel_t::PH>(g2_nn_cvt);
 template g2_nn_t identity<Channel_t::PH_bar>(g2_nn_cvt);
 template g2_nn_t identity<Channel_t::PP>(g2_nn_cvt);

 g2_nn_t identity_PH(g2_nn_vt g) { return identity<Channel_t::PH>(g); }
 g2_nn_t identity_PP(g2_nn_vt g) { return identity<Channel_t::PP>(g); }
 g2_nn_t identity_PH_bar(g2_nn_vt g) { return identity<Channel_t::PH_bar>(g); }
  
 template g2_iw_t identity<Channel_t::PH>(g2_iw_cvt);
 template g2_iw_t identity<Channel_t::PH_bar>(g2_iw_cvt);
 template g2_iw_t identity<Channel_t::PP>(g2_iw_cvt);

 g2_iw_t identity_PH(g2_iw_vt g) { return identity<Channel_t::PH>(g); }
 g2_iw_t identity_PP(g2_iw_vt g) { return identity<Channel_t::PP>(g); }
 g2_iw_t identity_PH_bar(g2_iw_vt g) { return identity<Channel_t::PH_bar>(g); }
  
} // namespace triqs_tprf
