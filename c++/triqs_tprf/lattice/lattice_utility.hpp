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

  g_w_t dlr_on_imfreq(g_Dc_cvt g_c, mesh::imfreq wmesh);
  g_t_t dlr_on_imtime(g_Dc_cvt g_c, mesh::imtime tmesh);
  chi_w_t dlr_on_imfreq(chi_Dc_cvt chi_c, mesh::imfreq wmesh);
  chi_t_t dlr_on_imtime(chi_Dc_cvt chi_c, mesh::imtime tmesh);

  /** Splits a rank 4 tensor-valued Green's function into dynamic and constant parts by tail fitting

    Splits a general rank 4 tensor-valued Green's function \f$ \chi_{abcd}(i\omega_n, \mathbf{k}) \f$ 
    into a dynamic and a constant part in Matsubara frequency space by fitting
    the high-frequency tail.
    
    \f[
        \chi_{abcd}(i\omega_n, \mathbf{k}) = 
            \chi^{(dyn)}_{abcd}(i\omega_n, \mathbf{k})
            + \chi^{(stat)}_{abcd}(\mathbf{k})
    \f]

  @param chi_wk : general rank 4 tensor-valued Green's function \f$ \chi_{abcd}(i\omega_n, \mathbf{k}) \f$. 
  @return Tuple of chi_dyn_wk, the dynamic part of chi \f$ \chi^{(dyn)}_{abcd}(i\omega_n, \mathbf{k}) \f$, which converges to zero for \f$ \omega_n \rightarrow \infty \f$, and chi_const_k, the part of chi that is constant in Matsubara frequency space \f$ \chi^{(stat)}_{abcd}(\mathbf{k}) \f$.
  */
  std::tuple<chi_wk_t, chi_k_t> split_into_dynamic_wk_and_constant_k(chi_wk_cvt chi_wk);

  /** Adds a dynamic and a static Green's function.

  @param g_dyn_fk : general Green's function \f$ G_{ab}(\omega, \mathbf{k}) \f$.
  @param g_stat_k : general Green's function \f$ G_{ab}(\mathbf{k}) \f$.
  @return g_fk : Green's function \f$ G_{ab}(\omega, \mathbf{k}) + G_{ab}(\mathbf{k}) \f$.
  */
  g_fk_t add_dynamic_and_static(g_fk_t g_dyn_fk, e_k_t g_stat_k);
  chi_fk_t add_dynamic_and_static(chi_fk_t chi_dyn_fk, chi_k_t chi_stat_k);

  /** Adds a dynamic and a static Green's function.

  @param g_dyn_wk : general Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$.
  @param g_stat_k : general Green's function \f$ G_{ab}(\mathbf{k}) \f$.
  @return g_wk : Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) + G_{ab}(\mathbf{k}) \f$.
  */
  g_wk_t add_dynamic_and_static(g_wk_t g_dyn_wk, e_k_t g_stat_k);
  chi_wk_t add_dynamic_and_static(chi_wk_t chi_dyn_wk, chi_k_t chi_stat_k);
  g_Dwk_t add_dynamic_and_static(g_Dwk_t g_dyn_wk, e_k_t g_stat_k);
  chi_Dwk_t add_dynamic_and_static(chi_Dwk_t chi_dyn_wk, chi_k_t chi_stat_k);

  /** Helper function to evaluate the Fermi-Dirac distribution function

  \f[
      f(\epsilon) = \frac{1}{\exp(\epsilon) + 1}
  \f]
  
  @param e : point at which to evaluate \f$ f(\epsilon) \f$.
  @return The value of \f$ f(\epsilon) \f$.
  */
  double fermi(double e);

  /** Helper function to evaluate the Bose-Einstein distribution function

  \f[
      n_B(\epsilon) = \frac{1}{\exp(\epsilon) - 1}
  \f]
  
  @param e : point at which to evaluate \f$ n_B(\epsilon) \f$.
  @return The value of \f$ n_B(\epsilon) \f$.
  */
  double bose(double e);
} // namespace triqs_tprf
