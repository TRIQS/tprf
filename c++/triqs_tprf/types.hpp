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
#pragma once

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <nda/nda.hpp>

namespace triqs_tprf {

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace nda;

/// Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)
enum class Channel_t { PP, PH, PH_bar };

typedef std::array<int,3> kpt_t;
  
/// Singl-particle response function type with one Matsubara frequency
typedef gf<imfreq, matrix_valued> g_iw_t;
typedef g_iw_t::const_view_type g_iw_cvt;
typedef g_iw_t::view_type g_iw_vt;

/// Lattice types
typedef gf<brzone, matrix_valued> ek_t;
typedef ek_t::const_view_type ek_cvt;
typedef ek_t::view_type ek_vt;
  
typedef gf<prod<imfreq, brzone>, matrix_valued> gk_iw_t;
typedef gk_iw_t::const_view_type gk_iw_cvt;
typedef gk_iw_t::view_type gk_iw_vt;

typedef gf<prod<imfreq, cyclat>, matrix_valued> gr_iw_t;
typedef gr_iw_t::const_view_type gr_iw_cvt;
typedef gr_iw_t::view_type gr_iw_vt;

typedef gf<prod<imtime, cyclat>, matrix_valued> gr_tau_t;
typedef gr_tau_t::const_view_type gr_tau_cvt;
typedef gr_tau_t::view_type gr_tau_vt;  

typedef gf<prod<imfreq, imfreq>, tensor_valued<4>> chi0_t;
typedef chi0_t::const_view_type chi0_cvt;
typedef chi0_t::view_type chi0_vt;

  // old style

typedef gf<prod<imfreq, imfreq, brzone>, tensor_valued<4>> chi0q_t;
typedef chi0q_t::const_view_type chi0q_cvt;
typedef chi0q_t::view_type chi0q_vt;
  
typedef gf<prod<imfreq, imfreq, cyclat>, tensor_valued<4>> chi0r_t;
typedef chi0r_t::const_view_type chi0r_cvt;
typedef chi0r_t::view_type chi0r_vt;

  // new style

typedef gf<prod<imfreq, imfreq, cyclat>, tensor_valued<4>> chi_wnr_t;
typedef chi_wnr_t::const_view_type chi_wnr_cvt;
typedef chi_wnr_t::view_type chi_wnr_vt;

typedef gf<prod<imfreq, imfreq, brzone>, tensor_valued<4>> chi_wnk_t;
typedef chi_wnk_t::const_view_type chi_wnk_cvt;
typedef chi_wnk_t::view_type chi_wnk_vt;

typedef gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> chi_wnn_t;
typedef chi_wnn_t::const_view_type chi_wnn_cvt;
typedef chi_wnn_t::view_type chi_wnn_vt;

typedef gf<prod<imfreq, imfreq>, tensor_valued<4>> chi_nn_t;
typedef chi_nn_t::const_view_type chi_nn_cvt;
typedef chi_nn_t::view_type chi_nn_vt;
  
typedef gf<prod<brzone, imfreq, imfreq, imfreq>, tensor_valued<4>> chi_kwnn_t;
typedef chi_kwnn_t::const_view_type chi_kwnn_cvt;
typedef chi_kwnn_t::view_type chi_kwnn_vt;
  
  // -- New style types FIXME notation
  
  // bosonic freq w for "\omega"
  // fermionic freq "n" for "\nu"
  // imagnary time "t" for "\tau"
  // real space "r"
  // momentum space "k"

typedef gf<brzone, matrix_valued> e_k_t;
typedef e_k_t::const_view_type e_k_cvt;
typedef e_k_t::view_type e_k_vt;

typedef gf<cyclat, matrix_valued> e_r_t;
typedef e_r_t::const_view_type e_r_cvt;
typedef e_r_t::view_type e_r_vt;

typedef gf<imfreq, matrix_valued> g_w_t;
typedef g_w_t::const_view_type g_w_cvt;
typedef g_w_t::view_type g_w_vt;

typedef gf<imtime, matrix_valued> g_t_t;
typedef g_t_t::const_view_type g_t_cvt;
typedef g_t_t::view_type g_t_vt;

typedef gf<prod<imfreq, brzone>, matrix_valued> g_wk_t;
typedef g_wk_t::const_view_type g_wk_cvt;
typedef g_wk_t::view_type g_wk_vt;

typedef gf<prod<imfreq, cyclat>, matrix_valued> g_wr_t;
typedef g_wr_t::const_view_type g_wr_cvt;
typedef g_wr_t::view_type g_wr_vt;

typedef gf<prod<imtime, cyclat>, matrix_valued> g_tr_t;
typedef g_tr_t::const_view_type g_tr_cvt;
typedef g_tr_t::view_type g_tr_vt;

typedef gf<prod<imtime, brzone>, matrix_valued> g_tk_t;
typedef g_tk_t::const_view_type g_tk_cvt;
typedef g_tk_t::view_type g_tk_vt;
  
typedef gf<imtime, tensor_valued<4>> chi_t_t;
typedef chi_t_t::const_view_type chi_t_cvt;
typedef chi_t_t::view_type chi_t_vt;

typedef gf<prod<imtime, cyclat>, tensor_valued<4>> chi_tr_t;
typedef chi_tr_t::const_view_type chi_tr_cvt;
typedef chi_tr_t::view_type chi_tr_vt;

typedef gf<prod<imfreq, cyclat>, tensor_valued<4>> chi_wr_t;
typedef chi_wr_t::const_view_type chi_wr_cvt;
typedef chi_wr_t::view_type chi_wr_vt;

typedef gf<prod<imfreq, brzone>, tensor_valued<4>> chi_wk_t;
typedef chi_wk_t::const_view_type chi_wk_cvt;
typedef chi_wk_t::view_type chi_wk_vt;

typedef gf<prod<brzone, imfreq>, tensor_valued<4>> chi_kw_t;
typedef chi_kw_t::const_view_type chi_kw_cvt;
typedef chi_kw_t::view_type chi_kw_vt;
  
typedef gf<brzone, tensor_valued<4>> chi_k_t;
typedef chi_k_t::const_view_type chi_k_cvt;
typedef chi_k_t::view_type chi_k_vt;

typedef gf<cyclat, tensor_valued<4>> chi_r_t;
typedef chi_r_t::const_view_type chi_r_cvt;
typedef chi_r_t::view_type chi_r_vt;

typedef gf<imfreq, tensor_valued<4>> chi_w_t;
typedef chi_w_t::const_view_type chi_w_cvt;
typedef chi_w_t::view_type chi_w_vt;

    // -- back to old style

typedef gf<prod<brzone, imfreq, imfreq, imfreq>, tensor_valued<4>> chiq_t;
typedef chiq_t::const_view_type chiq_cvt;
typedef chiq_t::view_type chiq_vt;

/// Two-particle response function type with three Matsubara frequencies
typedef gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_iw_t;
typedef g2_iw_t::const_view_type g2_iw_cvt;
typedef g2_iw_t::view_type g2_iw_vt;
  
// block greens functions
typedef block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> b_g2_iw_t;
typedef b_g2_iw_t::view_type b_g2_iw_vt;  
typedef b_g2_iw_t::const_view_type b_g2_iw_cvt;

typedef block_gf<imfreq, matrix_valued> b_g_iw_t;
typedef b_g_iw_t::view_type b_g_iw_vt;  
typedef b_g_iw_t::const_view_type b_g_iw_cvt;
  
/// Container type of one-particle Green and Vertex functions in Matsubara
//using g_iw_t = gf<imfreq, matrix_valued>;

/// Container type of two-particle Green and Vertex functions in imaginary time
using g_tau_t = gf<imtime, matrix_valued>;

/// Container type of $\chi_3$ in Matsubara frequencies
using chi2_iw_t = gf<imfreq, tensor_valued<4>>;

/// Container type of $\chi_3$ in imaginary time
using chi2_tau_t = gf<imtime, tensor_valued<4>>;

/// Container type of $\chi_3$ in Matsubara frequencies
using chi3_iw_t = gf<prod<imfreq, imfreq>, tensor_valued<4>>;

/// Container type of $\chi_3$ in imaginary time
using chi3_tau_t = gf<prod<imtime, imtime>, tensor_valued<4>>;

/// Container type of two-particle Green and Vertex functions in Matsubara
/// frequencies
using chi4_iw_t =
    gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>;

/// Container type of two-particle Green and Vertex functions in imaginary time
using chi4_tau_t =
    gf<prod<imtime, imtime, imtime>, tensor_valued<4>>;

} // namespace triqs_tprf
