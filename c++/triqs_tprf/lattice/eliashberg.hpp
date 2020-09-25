/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation and S. Käser
 * Authors: S. Käser, H. U.R. Strand
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

 /** Linearized Eliashberg product via summation

     Computes the linearized Eliashberg product in the singlet/triplet channel given by

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k}) 
         =
         -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
         \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
         \\
         \times
         G_{c\bar{e}}(i\nu',\mathbf{k}')
         G_{d\bar{f}}(-i\nu',-\mathbf{k}')
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{f}}(i\nu',\mathbf{k}')\,,

     by summation.

     @param Gamma_pp particle-particle vertex :math:`\Gamma^{\mathrm{s/t}}_{a\bar{b}c\bar{d}}(i\nu_n,\mathbf{k})`
     @param g_wk single particle Green's function :math:`G_{a\bar{b}}(i\nu_n,\mathbf{k})`
     @param delta_wk superconducting gap :math:`\Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})`
     @return Gives the result of the product :math:`\Delta^{\mathrm{s/t}, \mathrm{out}}`

  */

  g_wk_t eliashberg_product(chi_wk_vt Gamma_pp, g_wk_vt g_wk, g_wk_vt delta_wk);

 /** Linearized Eliashberg product via FFT

     Computes the linearized Eliashberg product in the singlet/triplet channel given by

     .. math::
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k}) 
        =
        -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
        \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
        \\
        \times
        G_{c\bar{e}}(i\nu',\mathbf{k}')
        G_{d\bar{f}}(-i\nu',-\mathbf{k}')
        \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{f}}(i\nu',\mathbf{k}')\,,

     by taking advantage of the convolution theorem.

     We therefore first calculate

     .. math::
         F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         G_{a\bar{c}}(i\nu,\mathbf{k})
         G_{b\bar{d}}(-i\nu,-\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})\,,

     which we then Fourier transform to imaginary time and real-space

     .. math::
        F^{\mathrm{s/t}}_{ab}(\tau,\mathbf{r})
        =
        \mathcal{F}^2
        \big(
        F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
        \big)\,.

     We then calculate first the dynamic gap
     
     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
         =
         -\frac{1}{2}
         \Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})
         F^{\mathrm{s/t}}_{cd}(\tau, \mathbf{r})\,,

     and then the static gap

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})
         =
         -\frac{1}{2}
         \Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})
         F^{\mathrm{s/t}}_{cd}(\tau=0, \mathbf{r})\,.

     We then Fourier transform the dynamic gap to imaginary frequencies

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        =
        \mathcal{F}
        \big(
        \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
        \big)\,,

     and then add both component together
     
     .. math::
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        =
        \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        +
        \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})\,,

    and then finally Fourier transform to :math:`\mathbf{k}`-space

    .. math::
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})
        =
        \mathcal{F}
        \big(
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        \big)\,.


     @param Gamma_pp_dyn_tr dynamic part of the particle-particle vertex :math:`\Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})`
     @param Gamma_pp_const_r static part of the particle-particle vertex :math:`\Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})`
     @param g_wk one-particle Green's function :math:`G_{a\bar{b}}(i\nu_n,\mathbf{k})`
     @param delta_wk superconducting gap :math:`\Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})`
     @return Gives the result of the product :math:`\Delta^{\mathrm{s/t}, \mathrm{out}}`

  */

  g_wk_t eliashberg_product_fft(chi_tr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt delta_wk);
  g_wk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt delta_wk);
  g_wk_t eliashberg_g_delta_g_product(g_wk_vt g_wk, g_wk_vt delta_wk);


  /** Split Gamma in dynamic and constant part by tail fitting

  @param Gamma_pp : particle-particle pairing vertex :math:`\Gamma(i\omega_n, \mathbf{k})`. 
  @return Tuple of Gamma_pp_dyn_wk, the dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`, and Gamma_pp_const_k, the part of Gamma that is constant in Matsubara frequency space :math:`\Gamma(\mathbf{k})`.
  */
  std::tuple<chi_wk_t, chi_k_t> split_into_dynamic_wk_and_constant_k(chi_wk_vt Gamma_pp);


  /** Fourier transform Gamma parts to imaginary time and real-space  
  
  @param Gamma_pp_dyn_wk : The dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`.
  @param Gamma_pp_const_k : The part of Gamma that is constant in Matsubara frequency space :math:`\Gamma(\mathbf{k})`.
  @return Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space, Gamma_pp_const_r, the constant part of Gamma in real-space.
  */
  std::tuple<chi_tr_t, chi_r_t> dynamic_and_constant_to_tr(chi_wk_vt Gamma_pp_dyn_wk, chi_k_vt Gamma_pp_const_k);
  e_r_t eliashberg_constant_gamma_f_product(chi_r_vt Gamma_pp_const_r, g_tr_t F_tr);

 /** The particle-particle vertex in the singlet channel

     Computes the singlet channel particle-particle vertex in the
     random phase approximation given by

     .. math::
         \Gamma^{\mathrm{singlet}}(i\omega_n,\mathbf{q}) =
         3 \mathbf{U}^{\mathrm{s}}
         \mathbf{\chi}^{\mathrm{s}}(i\omega_n,\mathbf{q})
         \mathbf{U}^{\mathrm{s}}
         -\mathbf{U}^{\mathrm{c}}
         \mathbf{\chi}^{\mathrm{c}}(i\omega_n,\mathbf{q})
         \mathbf{U}^{\mathrm{c}}
         + \frac{1}{2}\big(\mathbf{U}^{\mathrm{s}}+
         \mathbf{U}^{\mathrm{c}}\big)\,,

     where all products are particle-hole products.
     Note, that this is a special case, where the particle-particle vertex only
     depends on one bosonic frequency and momentum. It can therefore only be used
     in the linearized Eliashberg equation, if symmetries are enforced,
     as desribed in the theory here: :ref:`eliashberg_rpa`.

     @param chi_c charge susceptibility  :math:`\chi^{\mathrm{c}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q})`
     @param chi_s spin susceptibility  :math:`\chi^{\mathrm{s}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q})`
     @param U_c charge interaction  :math:`U^{\mathrm{c}}_{a\bar{b}c\bar{d}}`
     @param U_s spin interaction  :math:`U^{\mathrm{s}}_{a\bar{b}c\bar{d}}`
     @return The singlet channel particle-particle vertex :math:`\Gamma^{\mathrm{singlet}}(i\omega_n,\mathbf{q})`

  */

  chi_wk_t gamma_PP_singlet(chi_wk_vt chi_c, chi_wk_vt chi_s, array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s);

 /** The particle-particle vertex in the triplet channel

     Computes the triplet channel particle-particle vertex in the
     random phase approximation given by

     .. math::
         \Gamma^{\mathrm{triplet}}(i\omega_n,\mathbf{q}) =
         -\mathbf{U}^{\mathrm{s}}
         \mathbf{\chi}^{\mathrm{s}}(i\omega_n,\mathbf{q})
         \mathbf{U}^{\mathrm{s}}
         -\mathbf{U}^{\mathrm{c}}
         \mathbf{\chi}^{\mathrm{c}}(i\omega_n,\mathbf{q})
         \mathbf{U}^{\mathrm{c}}
         + \frac{1}{2}\big(\mathbf{U}^{\mathrm{s}}+
         \mathbf{U}^{\mathrm{c}}\big)\,,

     where all products are particle-hole products.
     Note, that this is a special case, where the particle-particle vertex only
     depends on one bosonic frequency and momentum. It can therefore only be used
     in the linearized Eliashberg equation, if symmetries are enforced,
     as desribed in the theory here: :ref:`eliashberg_rpa`.

     @param chi_c charge susceptibility  :math:`\chi^{\mathrm{c}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q})`
     @param chi_s spin susceptibility  :math:`\chi^{\mathrm{s}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q})`
     @param U_c charge interaction  :math:`U^{\mathrm{c}}_{a\bar{b}c\bar{d}}`
     @param U_s spin interaction  :math:`U^{\mathrm{s}}_{a\bar{b}c\bar{d}}`
     @return The triplet channel particle-particle vertex :math:`\Gamma^{\mathrm{triplet}}(i\omega_n,\mathbf{q})`

  */

  chi_wk_t gamma_PP_triplet(chi_wk_vt chi_c, chi_wk_vt chi_s, array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s);
  chi_wk_t gamma_PP_spin_charge(chi_wk_vt chi_c, chi_wk_vt chi_s, array_view<std::complex<double>, 4> U_c, array_view<std::complex<double>, 4> U_s, double charge_factor, double spin_factor);
}
