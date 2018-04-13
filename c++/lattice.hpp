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

#include "types.hpp"

#include <triqs/lattice/brillouin_zone.hpp>
#include <triqs/lattice/tight_binding.hpp>
#include <triqs/utility/mini_vector.hpp>

namespace tprf {

array<std::complex<double>, 6> cluster_mesh_fourier_interpolation(array<double, 2> k_vecs, chi_wr_cvt chi);

/** Construct a non-interacting Matsubara frequency lattice Green's function 

    Using a discretized dispersion $\epsilon(\mathbf{k})$, chemical potential $\mu$,
    and a Matsubara frequency Green's function mesh.

 .. math::
    G_{a\bar{b}}(k, i\omega_n) = \left[ 
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k} 
	\right]^{-1}_{a\bar{b}}

 @param mu chemical potential
 @param ek discretized lattice dispersion $\epsilon_{\bar{c}d}(\mathbf{k})$
 @param mesh imaginary frequency mesh
 */
gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh);
  
/** Construct an interacting Matsubara frequency lattice Green's function 

    Using a discretized dispersion $\epsilon_{\bar{a}b}(\mathbf{k})$, 
    chemical potential $\mu$, and a momentum independent Matsubara frequency 
    self energy $\Sigma_{\bar{a}b}(i\omega_n)$.

 .. math::
    G_{a\bar{b}}(k, i\omega_n) = \left[ 
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}}

 @param mu chemical potential
 @param ek discretized lattice dispersion $\epsilon_{\bar{c}d}(\mathbf{k})$
 @param sigma imaginary frequency self-energy $\Sigma_{\bar{a}b}(i\omega_n)$
 */
gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma);

gr_iw_t gr_from_gk(gk_iw_vt gk);
gk_iw_t gk_from_gr(gr_iw_vt gr);

gr_tau_t grt_from_grw(gr_iw_vt grw);

chi_tr_t chi0_tr_from_grt_PH(gr_tau_vt grt);
chi_wr_t chi_w0r_from_chi_tr(chi_tr_vt chi_tr);
chi_wr_t chi_wr_from_chi_tr(chi_tr_vt chi_tr, int nw);
chi_wk_t chi_wk_from_chi_wr(chi_wr_vt chi_wr);

chi_wk_t chi00_wk_from_ek(gf<brillouin_zone, matrix_valued> ek_in, int nw, double beta, double mu);
  
chi0r_t chi0r_from_gr_PH(int nw, int nnu, gr_iw_vt gr);

chi0r_t chi0r_from_chi0q(chi0q_vt chi0q);
chi0q_t chi0q_from_chi0r(chi0r_vt chi0r);

gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu(chi0q_t chi0q);
gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu_tail_corr_PH(chi0q_t chi0q);
gf<imfreq, tensor_valued<4>> chi0q_sum_nu_q(chi0q_t chi0q);

chiq_t chiq_from_chi0q_and_gamma_PH(chi0q_vt chi0q, g2_iw_vt gamma_ph);
gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu(chiq_t chiq);
gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(chiq_t chiq);
  
} // namespace tprf
