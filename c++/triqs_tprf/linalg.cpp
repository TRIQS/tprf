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

   auto g_inv = make_gf(g);

   using dat_t = array<g2_nn_cvt::scalar_t, g2_nn_cvt::data_rank, channel_memory_layout<CH>>;

   auto g_w_dat     = dat_t{g.data()};
   auto g_w_inv_dat = dat_t{g.data()};

   auto mat     = channel_matrix_view<CH>(g_w_dat);
   auto mat_inv = channel_matrix_view<CH>(g_w_inv_dat);

   mat_inv = inverse(mat);

   g_inv.data() = g_w_inv_dat;

   return g_inv;
 }
  
 /// Inverse: [G]^{-1}, Two-particle response-function inversion
 template <Channel_t CH> g2_iw_t inverse(g2_iw_cvt g) {

   //channel_grouping<CH> chg;
   auto g_inv = make_gf(g);

   for (auto const &w : std::get<0>(g.mesh())) {
     auto _         = all_t{};
     g_inv[w, _, _] = inverse<CH>(g[w, _, _]);
  }
  return g_inv;
 }

 // ----------------------------------------------------
  
 template <Channel_t CH> g2_nn_t product(g2_nn_cvt A, g2_nn_cvt B) {

   auto C = make_gf(A);

   using dat_t = array<g2_nn_cvt::scalar_t, g2_nn_cvt::data_rank, channel_memory_layout<CH>>;

   auto A_w_dat = dat_t{A.data()};
   auto B_w_dat = dat_t{B.data()};
   auto C_w_dat = dat_t{A_w_dat};

   auto A_mat = channel_matrix_view<CH>(A_w_dat);
   auto B_mat = channel_matrix_view<CH>(B_w_dat);
   auto C_mat = channel_matrix_view<CH>(C_w_dat);

   C_mat = A_mat * B_mat;

   C.data() = C_w_dat;

   return C;
 }

 /// product: C = A * B, two-particle response-function product
 template <Channel_t CH> g2_iw_t product(g2_iw_cvt A, g2_iw_cvt B) {

   auto C = make_gf(A);

   for (auto const &w : std::get<0>(A.mesh())) {
     auto _     = all_t{};
     C[w, _, _] = product<CH>(A[w, _, _], B[w, _, _]);
  }
  return C;
 }

 // ----------------------------------------------------

 template <Channel_t CH> g2_nn_t identity(g2_nn_cvt g) {

   using dat_t = array<g2_nn_cvt::scalar_t, g2_nn_cvt::data_rank, channel_memory_layout<CH>>;

   auto I = make_gf(g);

   auto I_w_dat = dat_t{I.data()};
   auto I_mat   = channel_matrix_view<CH>(I_w_dat);

   I_mat    = 1.0; // This sets a nda::matrix to the identity matrix...
   I.data() = I_w_dat;

   return I;
 }

 /// Identity: 1, identity two-particle response-function
 template <Channel_t CH> g2_iw_t identity(g2_iw_cvt g) {

   auto I = make_gf(g);

   for (auto const &w : std::get<0>(g.mesh())) {
     auto _     = all_t{};
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

 // ----------------------------------------------------

array<g2_nn_cvt::scalar_t, 4>
scalar_product_PH(g2_n_cvt vL, g2_nn_cvt M, g2_n_cvt vR) {

  using L_layout = contiguous_layout_with_stride_order<nda::encode(std::array<int, 5>({1, 2, 0, 4, 3}))>;
  using R_layout = contiguous_layout_with_stride_order<nda::encode(std::array<int, 5>({0, 1, 2, 3, 4}))>;
  using res_layout = contiguous_layout_with_stride_order<nda::encode(std::array<int, 4>({0, 1, 3, 2}))>;

  using L_vec_t = array<g2_n_cvt::scalar_t, g2_n_cvt::data_rank, L_layout>;
  using R_vec_t = array<g2_n_cvt::scalar_t, g2_n_cvt::data_rank, R_layout>;
  using mat_t = array<g2_nn_cvt::scalar_t, g2_nn_cvt::data_rank, channel_memory_layout<Channel_t::PH>>;
  
  array<g2_nn_cvt::scalar_t, 4, res_layout> res_dat(M.target_shape());

  auto M_dat = mat_t{M.data()};
  auto vL_dat = L_vec_t{vL.data()};
  auto vR_dat = R_vec_t{vR.data()};
  
  auto mat = channel_matrix_view<Channel_t::PH>(M_dat);
  matrix_view<g2_n_cvt::scalar_t> res(group_indices_view(res_dat, idx_group<0, 1>, idx_group<3, 2>));
  matrix_view<g2_n_cvt::scalar_t> vecL(group_indices_view(vL_dat, idx_group<1, 2>, idx_group<0, 4, 3>));
  matrix_view<g2_n_cvt::scalar_t> vecR(group_indices_view(vR_dat, idx_group<0, 1, 2>, idx_group<3, 4>)); 

  res = vecL * mat * vecR;

  return res_dat;
}

  
} // namespace triqs_tprf
