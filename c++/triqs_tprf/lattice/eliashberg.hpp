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

 /** Eliashberg product via summation

     Computes the Eliashberg product in the singlet/triplet channel given by

     \f[
         \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k}) 
         =
         -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
         \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
         F^{\mathrm{s/t}}_{dc}(i\nu',\mathbf{k}')\,,
     \f]

     by summation.

     Here \f$ F^{\mathrm{s/t}}_{ab}(i\nu_n,\mathbf{k}) \f$ is the anomalous propagator. 
     If linearized=true (default) it will evaluate the linearized propagator

     \f[
         F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         -G_{a\bar{c}}(i\nu,\mathbf{k})
         \bar{G}_{b\bar{d}}(i\nu,-\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})\,,
     \f]

     otherwise the full propagator will be evaluated (see :meth:`triqs_tprf.lattice.eliashberg_F_wk`).

     @param Gamma_pp particle-particle vertex \f$ \Gamma^{\mathrm{s/t}}_{a\bar{b}c\bar{d}}(i\nu_n,\mathbf{k}) \f$
     @param g_wk single particle Green's function \f$ G_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @param gbar_wk hole Green's function \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$ (optional; defaults to \f$ \bar{G}(i\nu_n,\mathbf{k})=-G(-i\nu_n,\mathbf{k}) \f$)
     @param delta_wk superconducting gap \f$ \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @param linearized: a boolean determining whether the linearized Eliashberg product should be returned (optional; defaults to true)
     @return Gives the result of the product \f$ \Delta^{\mathrm{s/t}, \mathrm{out}} \f$

  */

  g_wk_t eliashberg_product(chi_wk_vt Gamma_pp, g_wk_vt g_wk, g_wk_vt delta_wk, bool linearized=true);
  g_wk_t eliashberg_product(chi_wk_vt Gamma_pp, g_wk_vt g_wk, g_wk_vt gbar_wk, g_wk_vt delta_wk, bool linearized=true);

 /** Eliashberg product via FFT

     Computes the Eliashberg product in the singlet/triplet channel given by

     \f[
         \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k}) 
         =
         -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
         \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
         \\
         \times
         F^{\mathrm{s/t}}_{dc}(i\nu',\mathbf{k}')\,,
     \f]

     by taking advantage of the convolution theorem.

     We therefore first calculate the anomalous propagator \f$ F^{\mathrm{s/t}}_{ab}(i\nu_n,\mathbf{k}) \f$.
     If linearized=true (default) we will evaluate the linearized propagator

     \f[
         F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         -G_{a\bar{c}}(i\nu,\mathbf{k})
         \bar{G}_{b\bar{d}}(i\nu,-\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})\,,
     \f]

     otherwise the full propagator will be evaluated (see :meth:`triqs_tprf.lattice.eliashberg_F_wk`).

     We then Fourier transform it to imaginary time and real-space

     \f[
        F^{\mathrm{s/t}}_{ab}(\tau,\mathbf{r})
        =
        \mathcal{F}^2
        \big(
        F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
        \big)\,.
     \f]

     We then calculate first the dynamic gap
     
     \f[
         \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
         =
         -\frac{1}{2}
         \Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})
         F^{\mathrm{s/t}}_{cd}(\tau, \mathbf{r})\,,
     \f]

     and then the static gap

     \f[
         \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})
         =
         -\frac{1}{2}
         \Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})
         F^{\mathrm{s/t}}_{cd}(\tau=0, \mathbf{r})\,.
     \f]

     We then Fourier transform the dynamic gap to imaginary frequencies

     \f[
         \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        =
        \mathcal{F}
        \big(
        \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
        \big)\,,
     \f]

     and then add both component together
     
     \f[
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        =
        \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        +
        \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})\,,
     \f]

    and then finally Fourier transform to \f$ \mathbf{k} \f$-space

    \f[
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})
        =
        \mathcal{F}
        \big(
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        \big)\,.
    \f]


     @param Gamma_pp_dyn_tr dynamic part of the particle-particle vertex \f$ \Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r}) \f$
     @param Gamma_pp_const_r static part of the particle-particle vertex \f$ \Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r}) \f$
     @param g_wk one-particle Green's function \f$ G_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @param gbar_wk hole Green's function \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$ (optional; defaults to \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k})=-G_{a\bar{b}}(-i\nu_n,\mathbf{k}) \f$)
     @param delta_wk superconducting gap \f$ \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @param linearized: a boolean determining whether the linearized Eliashberg product should be returned (optional; defaults to true)
     @return Gives the result of the product \f$ \Delta^{\mathrm{s/t}, \mathrm{out}} \f$

  */

  g_wk_t eliashberg_product_fft(chi_tr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt delta_wk, bool linearized=true);
  g_wk_t eliashberg_product_fft(chi_tr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt gbar_wk, g_wk_vt delta_wk, bool linearized=true);
  g_Dwk_t eliashberg_product_fft(chi_Dtr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_Dwk_vt g_wk, g_Dwk_vt delta_wk, bool linearized=true);
  g_Dwk_t eliashberg_product_fft(chi_Dtr_vt Gamma_pp_dyn_tr, chi_r_vt Gamma_pp_const_r, g_Dwk_vt g_wk, g_Dwk_vt gbar_wk, g_Dwk_vt delta_wk, bool linearized=true);

  g_wk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt gbar_wk, g_wk_vt delta_wk, bool linearized=true);
  g_wk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_wk_vt g_wk, g_wk_vt delta_wk, bool linearized=true);
  g_Dwk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_Dwk_vt g_wk, g_Dwk_vt gbar_wk, g_Dwk_vt delta_wk, bool linearized=true);
  g_Dwk_t eliashberg_product_fft_constant(chi_r_vt Gamma_pp_const_r, g_Dwk_vt g_wk, g_Dwk_vt delta_wk, bool linearized=true);
  g_wk_t eliashberg_g_delta_g_product(g_wk_vt g_wk, g_wk_vt gbar_wk, g_wk_vt delta_wk);
  g_wk_t eliashberg_g_delta_g_product(g_wk_vt g_wk, g_wk_vt delta_wk);
  g_Dwk_t eliashberg_g_delta_g_product(g_Dwk_vt g_wk, g_Dwk_vt gbar_wk, g_Dwk_vt delta_wk);
  g_Dwk_t eliashberg_g_delta_g_product(g_Dwk_vt g_wk, g_Dwk_vt delta_wk);

 /** Anomalous propagator in Eliashberg theory

     Computes the full (i.e., not linearized) anomalous propagator in Eliashberg theory.
     We first construct the denominator

     \f[
         \Theta^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         \delta_{ab}
         -
         G_{a\bar{c}}(i\nu,\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})
         \bar{G}_{e\bar{d}}(i\nu,-\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{b}}(i\nu,\mathbf{k})\,,
     \f]

     which we then invert and substitute into

     \f[
         F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         -\Theta^{\mathrm{s/t}, -1}_{ac}(i\nu,\mathbf{k})
         G_{c\bar{d}}(i\nu,\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{d}\bar{e}}(i\nu,\mathbf{k})
         \bar{G}_{b\bar{e}}(-i\nu,-\mathbf{k})\,.
     \f]

     WARNING: This method has only been tested for single-orbital systems.
     For multi-orbital systems the index ordering should be debugged carefully.

     @param g_wk single particle Green's function \f$ G_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @param gbar_wk hole Green's function \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$ (optional; defaults to \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k})=-G_{a\bar{b}}(-i\nu_n,\mathbf{k}) \f$)
     @param delta_wk superconducting gap \f$ \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @return The anomalous propagator \f$ F^{\mathrm{s/t}}_{ab}(i\nu_n,\mathbf{k}) \f$

  */

  g_wk_t eliashberg_F_wk(g_wk_vt g_wk, g_wk_vt gbar_wk, g_wk_vt delta_wk);
  g_wk_t eliashberg_F_wk(g_wk_vt g_wk, g_wk_vt delta_wk);
  g_Dwk_t eliashberg_F_wk(g_Dwk_vt g_wk, g_Dwk_vt gbar_wk, g_Dwk_vt delta_wk);
  g_Dwk_t eliashberg_F_wk(g_Dwk_vt g_wk, g_Dwk_vt delta_wk);


 /** Normal propagator in Eliashberg theory

     Computes the full (i.e., not linearized) normal propagator in Eliashberg theory.
     We first construct the denominator

     \f[
         \Theta^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         \delta_{ab}
         -
         G_{a\bar{c}}(i\nu,\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})
         \bar{G}_{e\bar{d}}(i\nu,-\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{b}}(i\nu,\mathbf{k})\,,
     \f]

     which we then invert and substitute into

     \f[
         G^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         \Theta^{\mathrm{s/t}, -1}_{ac}(i\nu,\mathbf{k})
         G_{c\bar{b}}(i\nu,\mathbf{k})
     \f]

     WARNING: This method has only been tested for single-orbital systems.
     For multi-orbital systems the index ordering should be debugged carefully.

     @param g_in_wk single particle Green's function \f$ G_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @param gbar_wk hole Green's function \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k}) \f$ (optional; defaults to \f$ \bar{G}_{a\bar{b}}(i\nu_n,\mathbf{k})=-G_{a\bar{b}}(-i\nu_n,\mathbf{k}) \f$)
     @param delta_wk superconducting gap \f$ \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k}) \f$
     @return The normal propagator \f$ G^{\mathrm{s/t}}_{ab}(i\nu_n,\mathbf{k}) \f$

  */

  g_wk_t eliashberg_g_wk(g_wk_vt g_in_wk, g_wk_vt gbar_wk, g_wk_vt delta_wk);
  g_wk_t eliashberg_g_wk(g_wk_vt g_in_wk, g_wk_vt delta_wk);
  g_Dwk_t eliashberg_g_wk(g_Dwk_vt g_in_wk, g_Dwk_vt gbar_in_wk, g_Dwk_vt delta_wk);
  g_Dwk_t eliashberg_g_wk(g_Dwk_vt g_in_wk, g_Dwk_vt delta_wk);


  /** Fourier transform Gamma parts to imaginary time and real-space  
  
  @param Gamma_pp_dyn_wk : The dynamic part of Gamma, which converges to zero for \f$ \omega_n \rightarrow \infty \f$.
  @param Gamma_pp_const_k : The part of Gamma that is constant in Matsubara frequency space \f$ \Gamma(\mathbf{k}) \f$.
  @return Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for \f$ \omega_n \rightarrow \infty \f$, but now in \f$ \tau \f$-space, Gamma_pp_const_r, the constant part of Gamma in real-space.
  */
  std::tuple<chi_tr_t, chi_r_t> dynamic_and_constant_to_tr(chi_wk_vt Gamma_pp_dyn_wk, chi_k_vt Gamma_pp_const_k);

  /** Fourier transform Gamma parts to imaginary time and real-space  
  
  @param Gamma_pp_dyn_wk : The dynamic part of Gamma, which converges to zero for \f$ \omega_n \rightarrow \infty \f$.
  @param Gamma_pp_const_k : The part of Gamma that is constant in Matsubara frequency space \f$ \Gamma(\mathbf{k}) \f$.
  @return Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for \f$ \omega_n \rightarrow \infty \f$, but now in \f$ \tau \f$-space, Gamma_pp_const_r, the constant part of Gamma in real-space.
  */
  std::tuple<chi_Dtr_t, chi_r_t> dynamic_and_constant_to_tr(chi_Dwk_vt Gamma_pp_dyn_wk, chi_k_vt Gamma_pp_const_k);

  e_r_t eliashberg_constant_gamma_f_product(chi_r_vt Gamma_pp_const_r, g_tr_t F_tr);


  /** Computes reducible ladder vertex for the approximation of a local and static vertex.

    In this approximation the reducible ladder vertex in density/magnetic channel are given by

    \f[
        \Phi^{\text{d/m}}_{a\overline{b}c\overline{d}}(Q)
        &\approx
        \frac{1}{(N_\mathbf{k}\beta)^2}
        \sum_{K'', K'''}
        U^{\text{d/m}}\chi^{\text{d/m}}(Q, K'', K''') U^{\text{d/m}}
        \\
        &\approx
        U^{\mathrm{d/m}}
        \chi^{\text{d/m}}(Q) U^{\mathrm{d/m}}\,,
    \f]


    where all products are particle-hole products.
    The reducible ladder vertex in then only dependent on one bosonic frequency and momentum.
    It can then be used in :meth:`triqs_tprf.eliashberg.construct_gamma_singlet_rpa`
    or :meth:`triqs_tprf.eliashberg.construct_gamma_rpa` to construct the
    irreducible singlet/triplet vertex.

    @param chi density/magnetic susceptibility  \f$ \chi^{\mathrm{d/m}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q}) \f$
    @param U density/magnetic local and static vertex  \f$ U^{\mathrm{d/m}}_{a\bar{b}c\bar{d}} \f$
    @return The reducible ladder vertex in the density/magnetic channel \f$ \Phi^{\mathrm{d/m}}(i\omega_n,\mathbf{q}) \f$

  */
  chi_wk_t construct_phi_wk(chi_wk_vt chi, array_contiguous_view<std::complex<double>, 4> U);
}
