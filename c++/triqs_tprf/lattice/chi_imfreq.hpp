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

#include "../types.hpp"

namespace triqs_tprf {

  /** Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})`.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})

  @param nw Number of bosonic Matsubara freqiencies.
  @param nn Number of fermionic Matsubara freqiencies.
  @param g_tr Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.
 */
  chi_wnr_t chi0r_from_gr_PH(int nw, int nn, g_wr_cvt g_nr);

  /** Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` without MPI parallellization.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})

  @param nw Number of bosonic Matsubara freqiencies.
  @param nn Number of fermionic Matsubara freqiencies.
  @param g_tr Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.
 */
  chi_wnr_t chi0r_from_gr_PH_nompi(int nw, int nn, g_wr_cvt g_nr);

  /** Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` with convolution in k-space.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
     - \frac{\beta}{N_k} \sum_\mathbf{k} 
     G_{d\bar{a}}(\nu, \mathbf{k}) \cdot G_{b\bar{c}}(\nu + \omega, \mathbf{k} - \mathbf{q})

  @param nw Number of bosonic Matsubara freqiencies.
  @param nn Number of fermionic Matsubara freqiencies.
  @param g_tr Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.
 */
  chi_wnk_t chi0q_from_g_wk_PH(int nw, int nn, g_wk_cvt g_wk);

  /** Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum-space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real-space.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     \mathcal{F}_{\mathbf{q} \rightarrow \mathbf{r}} \left\{
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})
     \right\}

  @param chi_wnk Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.
 */
  chi_wnr_t chi0r_from_chi0q(chi_wnk_cvt chi_wnk);

  /** Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum space.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
     \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{q}} \left\{
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})
     \right\}

  @param chi_wnr Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.
  @return Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.
 */
  chi_wnk_t chi0q_from_chi0r(chi_wnr_cvt chi_wnr);

  /** Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max} \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

  @param chi_wnk Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.
  @return Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.
 */
  chi_wk_t chi0q_sum_nu(chi_wnk_cvt chi_wnk);

  /** Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})` using higher order tail corrections when summing to infinity.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{\beta^2} \sum_{\nu=-\infty}^\infty \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

  @param chi_wnk Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.
  @return Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.
 */
  chi_wk_t chi0q_sum_nu_tail_corr_PH(chi_wnk_cvt chi_wnk);

  /** Sum over fermionic frequency and momentum in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{N_k} \sum_\matbf{k} \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max}
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

  @param chi_wnk Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.
  @return Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega)` in one bosonic Matsubara frequency.
 */
  chi_w_t chi0q_sum_nu_q(chi_wnk_cvt chi_wnk);

  /** Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k})`.

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k}) =
     \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}

  @param chi0_wnk Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
  @param gamma_ph_wnn Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.
  @return Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k})`.
 */
  chi_kwnn_t chiq_from_chi0q_and_gamma_PH(chi_wnk_cvt chi0_wnk, chi_wnn_cvt gamma_ph_wnn);

  /** Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}

  @param chi0_wnk Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
  @param gamma_ph_wnn Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.
  @return Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
 */
  chi_kw_t chiq_sum_nu_from_chi0q_and_gamma_PH(chi_wnk_cvt chi0_wnk, chi_wnn_cvt gamma_ph_wnn);

  chi_kwnn_t f_from_chi0q_and_gamma_PH(chi_wnk_cvt chi0_wnk, chi_wnn_cvt gamma_ph_wnn);

  gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu_from_g_wk_and_gamma_PH(gk_iw_t g_wk, g2_iw_vt gamma_ph_wnn,
                                                                                                     int tail_corr_nwf = -1);

  gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>>
  chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH(double mu, ek_vt e_k, g_iw_vt sigma_w, g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1);

  gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu(chiq_t chiq);

  gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(chiq_t chiq);

} // namespace triqs_tprf
