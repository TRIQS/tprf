/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation and S. Käser
 * Authors: H. U.R. Strand, S. Käser
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

namespace tprf {

 /** Linearized Eliashberg product

     Computes the product

     .. math::
         \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{k}'} \sum_{i\nu'}
	 \Gamma_{A\bar{a}B\bar{b}}(\mathbf{k}-\mathbf{k}', i\nu - i\nu')
	 \\ \times
	 G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')

     @param chi_pp particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\nu_n)`
     @param g_kw single particle Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\nu_n)`
     @param delta_kw pairing self-energy :math:`\Delta_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n)`
     @return Gives the result of the product :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG \Delta`

  */

  gk_iw_t eliashberg_product(chi_wk_vt Gamma_pp, gk_iw_vt g_wk, gk_iw_vt delta_wk);

 /** Linearized Eliashberg product via FFT

     Computes the product

     .. math::
         \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{k}'} \sum_{i\nu'}
	 \Gamma_{A\bar{a}B\bar{b}}(\mathbf{k}-\mathbf{k}', i\nu - i\nu')
	 \\ \times
	 G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')\,,

     by taking advantage of the convolution theorem.

     We therefore first calculate

     .. math::
        \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{r}, \tau) = 
	 -\Gamma_{A\bar{a}B\bar{b}}(\mathbf{r}, \tau) F_{AB}(\mathbf{r}, \tau) \,,

     where 

     .. math::
        F_{AB}(\mathbf{r}, \tau)  = 
        \mathcal{F}\big(G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')\big)\,.

     Then we Fourier transform 

     .. math::
          \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) = 
          \mathcal{F}\big(\Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{r}, \tau)\big)\,,

    to get the same result, but with far less computational effort.

     @param chi_rt dynamic part of the particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{r}, \tau)`
     @param chi_r constant part of the particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{r})`
     @param g_kw single particle Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\nu_n)`
     @param delta_kw pairing self-energy :math:`\Delta_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n)`
     @return Gives the result of the product :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG \Delta`

  */

  gk_iw_t eliashberg_product_fft(chi_tr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, gk_iw_vt g_wk, gk_iw_vt delta_wk);
  gk_iw_t eliashberg_g_delta_g_product(gk_iw_vt g_wk, gk_iw_vt delta_wk);
  std::tuple<chi_wk_vt, chi_k_vt> split_into_dynamic_wk_and_constant_k(chi_wk_vt Gamma_pp);
  std::tuple<chi_tr_vt, chi_r_vt> dynamic_and_constant_to_tr(chi_wk_vt Gamma_pp_dyn_wk, chi_k_vt Gamma_pp_const_k);

 /** Gamma particle-particle singlet

     Computes the particle-particle vertex for singlet pairing in the RPA limit

    .. math::
        \Gamma^{(\mathrm{singlet})}(a\bar{b}c\bar{d}) =
        \frac{3}{2} U^{(\mathrm{s})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{s})}(\bar{B}A\bar{C}D) 
        U^{(\mathrm{s})}(D\bar{C}c\bar{d}) \\
        -\frac{1}{2} U^{(\mathrm{c})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{c})}(\bar{B}A\bar{C}D) 
        U^{(\mathrm{c})}(D\bar{C}c\bar{d}) \\
       + \frac{1}{2}\big(U^{(\mathrm{s})}(a\bar{b}c\bar{d})+
        U^{(\mathrm{c})}(a\bar{b}c\bar{d})\big)

     @param chi_c charge susceptibility  :math:`\chi^{(\mathrm{c})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`
     @param chi_s spin susceptibility  :math:`\chi^{(\mathrm{s})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`
     @param U_c charge interaction  :math:`U^{(\mathrm{c})}_{a\bar{b}c\bar{d}}`
     @param U_s spin interaction  :math:`U^{(\mathrm{s})}_{a\bar{b}c\bar{d}}`
     @return :math:`\Gamma^{(\mathrm{singlet})}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\omega_n)`

  */

  chi_wk_t gamma_PP_singlet(chi_wk_vt chi_c, chi_wk_vt chi_s, array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s);

 /** Gamma particle-particle triplet

     Computes the particle-particle vertex for triplet pairing in the RPA limit

    .. math::
        \Gamma^{(\mathrm{triplet})}(a\bar{b}c\bar{d}) =
        -\frac{1}{2} U^{(\mathrm{s})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{s})}(\bar{B}A\bar{C}D) 
        U^{(\mathrm{s})}(D\bar{C}c\bar{d}) \\
        -\frac{1}{2} U^{(\mathrm{c})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{c})}(\bar{B}A\bar{C}D) 
        U^{(\mathrm{c})}(D\bar{C}c\bar{d}) \\
       + \frac{1}{2}\big(U^{(\mathrm{s})}(a\bar{b}c\bar{d})+
        U^{(\mathrm{c})}(a\bar{b}c\bar{d})\big)

     @param chi_c charge susceptibility  :math:`\chi^{(\mathrm{c})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`
     @param chi_s spin susceptibility  :math:`\chi^{(\mathrm{s})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`
     @param U_c charge interaction  :math:`U^{(\mathrm{c})}_{a\bar{b}c\bar{d}}`
     @param U_s spin interaction  :math:`U^{(\mathrm{s})}_{a\bar{b}c\bar{d}}`
     @return :math:`\Gamma^{(\mathrm{triplet})}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\omega_n)`

  */

  chi_wk_t gamma_PP_triplet(chi_wk_vt chi_c, chi_wk_vt chi_s, array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s);
}
