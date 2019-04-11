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
#include <triqs/arrays.hpp>

namespace tprf {

using namespace triqs::gfs;
using namespace triqs::arrays;

/// Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)
enum class Channel_t { PP, PH, PH_bar };

typedef triqs::utility::mini_vector<int,3> kpt_t;
  
/// Singl-particle response function type with one Matsubara frequency
typedef gf<imfreq, matrix_valued> g_iw_t;
typedef g_iw_t::const_view_type g_iw_cvt;
typedef g_iw_t::view_type g_iw_vt;

/// Lattice types
typedef gf<brillouin_zone, matrix_valued> ek_t;
typedef ek_t::const_view_type ek_cvt;
typedef ek_t::view_type ek_vt;
  
typedef gf<cartesian_product<imfreq, brillouin_zone>, matrix_valued> gk_iw_t;
typedef gk_iw_t::const_view_type gk_iw_cvt;
typedef gk_iw_t::view_type gk_iw_vt;

typedef gf<cartesian_product<imfreq, cyclic_lattice>, matrix_valued> gr_iw_t;
typedef gr_iw_t::const_view_type gr_iw_cvt;
typedef gr_iw_t::view_type gr_iw_vt;

typedef gf<cartesian_product<imtime, cyclic_lattice>, matrix_valued> gr_tau_t;
typedef gr_tau_t::const_view_type gr_tau_cvt;
typedef gr_tau_t::view_type gr_tau_vt;  

typedef gf<cartesian_product<imfreq, imfreq>, tensor_valued<4>> chi0_t;
typedef chi0_t::const_view_type chi0_cvt;
typedef chi0_t::view_type chi0_vt;

typedef gf<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q_t;
typedef chi0q_t::const_view_type chi0q_cvt;
typedef chi0q_t::view_type chi0q_vt;

typedef gf<cartesian_product<imfreq, imfreq, cyclic_lattice>, tensor_valued<4>> chi0r_t;
typedef chi0r_t::const_view_type chi0r_cvt;
typedef chi0r_t::view_type chi0r_vt;

  // -- New style types FIXME notation
  
  // bosonic freq w for "\omega"
  // fermionic freq "n" for "\nu"
  // imagnary time "t" for "\tau"
  // real space "r"
  // momentum space "k"

typedef gf<cartesian_product<imfreq, brillouin_zone>, matrix_valued> g_wk_t;
typedef gk_iw_t::const_view_type g_wk_cvt;
typedef gk_iw_t::view_type g_wk_vt;
  
typedef gf<imtime, tensor_valued<4>> chi_t_t;
typedef chi_t_t::const_view_type chi_t_cvt;
typedef chi_t_t::view_type chi_t_vt;

typedef gf<cartesian_product<imtime, cyclic_lattice>, tensor_valued<4>> chi_tr_t;
typedef chi_tr_t::const_view_type chi_tr_cvt;
typedef chi_tr_t::view_type chi_tr_vt;

typedef gf<cartesian_product<imfreq, cyclic_lattice>, tensor_valued<4>> chi_wr_t;
typedef chi_wr_t::const_view_type chi_wr_cvt;
typedef chi_wr_t::view_type chi_wr_vt;

typedef gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi_wk_t;
typedef chi_wk_t::const_view_type chi_wk_cvt;
typedef chi_wk_t::view_type chi_wk_vt;

typedef gf<brillouin_zone, tensor_valued<4>> chi_k_t;
typedef chi_k_t::const_view_type chi_k_cvt;
typedef chi_k_t::view_type chi_k_vt;

typedef gf<cyclic_lattice, tensor_valued<4>> chi_r_t;
typedef chi_r_t::const_view_type chi_r_cvt;
typedef chi_r_t::view_type chi_r_vt;

    // -- back to old style

typedef gf<cartesian_product<brillouin_zone, imfreq, imfreq, imfreq>, tensor_valued<4>> chiq_t;
typedef chiq_t::const_view_type chiq_cvt;
typedef chiq_t::view_type chiq_vt;

/// Two-particle response function type with three Matsubara frequencies
typedef gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_iw_t;

/// Short hand constant view type for g2_iw_t
typedef g2_iw_t::const_view_type g2_iw_cvt;

/// Short hand view type for g2_iw_t
typedef g2_iw_t::view_type g2_iw_vt;
  
using scalar_t = g2_iw_t::scalar_t;

// block greens functions
typedef block2_gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> b_g2_iw_t;
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
using chi3_iw_t = gf<cartesian_product<imfreq, imfreq>, tensor_valued<4>>;

/// Container type of $\chi_3$ in imaginary time
using chi3_tau_t = gf<cartesian_product<imtime, imtime>, tensor_valued<4>>;

/// Container type of two-particle Green and Vertex functions in Matsubara
/// frequencies
using chi4_iw_t =
    gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;

/// Container type of two-particle Green and Vertex functions in imaginary time
using chi4_tau_t =
    gf<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>;

} // namespace tprf
