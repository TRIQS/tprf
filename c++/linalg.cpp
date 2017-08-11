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

namespace tprf {

 // ----------------------------------------------------

 /// Inverse: [G]^{-1}, Two-particle response-function inversion
 template <Channel_t CH> g2_iw_t inverse(g2_iw_t const &g) {

  channel_grouping<CH> chg;
  auto g_inv = make_gf(g.mesh(), g.target(), g.memory_layout());

  for (auto const &w : std::get<0>(g.mesh())) {
   auto _ = var_t{};

   auto g_w = make_gf(g[w][_][_], chg.memory_layout());
   auto g_w_inv = make_gf(g_w.mesh(), g_w.target(), g_w.memory_layout());

   auto mat = chg.matrix_view(g_w.data());
   auto mat_inv = chg.matrix_view(g_w_inv.data());

   mat_inv = inverse(mat);

   g_inv[w][_][_] = g_w_inv;
  }
  return g_inv;
 }

 // ----------------------------------------------------

 /// product: C = A * B, two-particle response-function product
 template <Channel_t CH> g2_iw_t product(g2_iw_t const &A, g2_iw_t const &B) {

  channel_grouping<CH> chg;
  // check that A and B are compatible!

  auto C = make_gf(A.mesh(), A.target(), A.memory_layout());

  for (auto const &w : std::get<0>(A.mesh())) {
   auto _ = var_t{};

   auto A_w = make_gf(A[w][_][_], chg.memory_layout());
   auto B_w = make_gf(B[w][_][_], chg.memory_layout());
   auto C_w = make_gf(A_w.mesh(), A_w.target(), A_w.memory_layout());

   auto A_mat = chg.matrix_view(A_w.data());
   auto B_mat = chg.matrix_view(B_w.data());
   auto C_mat = chg.matrix_view(C_w.data());

   C_mat = A_mat * B_mat;

   C[w][_][_] = C_w;
  }
  return C;
 }

 // ----------------------------------------------------

 /// Identity: 1, identity two-particle response-function
 template <Channel_t CH> g2_iw_t identity(g2_iw_t const &g) {

  channel_grouping<CH> chg;
  auto I = make_gf(g.mesh(), g.target(), g.memory_layout());

  for (auto const &w : std::get<0>(g.mesh())) {
   auto _ = var_t{};
   auto I_w = make_gf(I[w][_][_], chg.memory_layout());
   auto I_mat = chg.matrix_view(I_w.data());

   I_mat = 1.0; // This sets a triqs::arrays::matrix to the identity matrix...
   I[w][_][_] = I_w;
  }
  return I;
 }

 // ----------------------------------------------------

 template g2_iw_t inverse<Channel_t::PH>(g2_iw_t const &);
 template g2_iw_t product<Channel_t::PH>(g2_iw_t const &, g2_iw_t const &);
 template g2_iw_t identity<Channel_t::PH>(g2_iw_t const &);

} // namespace tprf
