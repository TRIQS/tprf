/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, The Simons Foundation
 * Authors: H. U.R. Strand, Y. in 't Veld, M. Rösner
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

  /** Density matrix from lattic Green's function
      
    @param g_wk single particle Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$
    @return rho_k density matrix \f$ \rho_{ab}(\mathbf{k}) \f$
  */
  e_k_t rho_k_from_g_wk(g_wk_cvt g_wk);
  e_k_t rho_k_from_g_wk(g_Dwk_cvt g_wk);

  /** GW self energy \f$ \Sigma(i\omega_n, \mathbf{k}) \f$ calculator for dynamic interactions

    Splits the interaction into a dynamic and a static part
    
    \f[
        W_{abcd}(i\omega_n, \mathbf{k}) = 
            W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k})
            + V_{abcd}(\mathbf{k})
    \f]

    by fitting the high-frequency tail.

    Fourier transforms the dynamic part of the interaction and the 
    single-particle Green's function to imaginary time and real space.

    \f[
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}
    \f]

    \f[
        W^{(dyn)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k}) \right\}
    \f]

    computes the GW self-energy as the product

    \f[
        \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) =
          - \sum_{cd} W^{(dyn)}_{acdb}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})
    \f]

    and transforms back to frequency and momentum

    \f[
        \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) \right\}
    \f]

    The self-energy of the static part of the interaction is calculated
    as the sum

    \f[
        \Sigma^{(stat)}_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{k}) \rho_{dc}(\mathbf{k} + \mathbf{q})
    \f]

    where \f$ \rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k}) \f$ is the density matrix of the
    single particle Green's function.

    The total GW self-energy is given by

    \f[
        \Sigma_{ab}(i\omega_n, \mathbf{k}) = 
          \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k})
          + \Sigma^{(stat)}_{ab}(\mathbf{k})
    \f]

    @param W_wk interaction \f$ W_{abcd}(i\omega_n, \mathbf{k}) \f$
    @param g_wk single particle Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$
    @return GW self-energy \f$ \Sigma_{ab}(i\omega_n, \mathbf{k}) \f$
 */

  g_wk_t gw_sigma(chi_wk_cvt W_wk, g_wk_cvt g_wk);

  /** GW self energy \f$ \Sigma(i\omega_n, \mathbf{k}) \f$ calculator for dynamic interactions

    Fourier transforms the dynamic part of the interaction and the 
    single-particle Green's function to imaginary time and real space.

    \f[
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}
    \f]

    \f[
        W^{(dyn)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k}) \right\}
    \f]

    computes the GW self-energy as the product

    \f[
        \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) =
          - \sum_{cd} W^{(dyn)}_{acdb}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})
    \f]

    and transforms back to frequency and momentum

    \f[
        \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) \right\}
    \f]

    The self-energy of the static part of the interaction is calculated
    as the sum

    \f[
        \Sigma^{(stat)}_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{k}) \rho_{dc}(\mathbf{k} + \mathbf{q})
    \f]

    where \f$ \rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k}) \f$ is the density matrix of the
    single particle Green's function.

    The total GW self-energy is given by

    \f[
        \Sigma_{ab}(i\omega_n, \mathbf{k}) = 
          \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k})
          + \Sigma^{(stat)}_{ab}(\mathbf{k})
    \f]

    @param W_wk interaction \f$ W_{abcd}(i\omega_n, \mathbf{k}) \f$
    @param v_k static interaction \f$ V_{abcd}(\mathbf{q}) \f$
    @param g_wk single particle Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$
    @return GW self-energy \f$ \Sigma_{ab}(i\omega_n, \mathbf{k}) \f$
 */

  g_Dwk_t gw_sigma(chi_Dwk_cvt W_wk, chi_k_cvt v_k, g_Dwk_cvt g_wk);

  /** Hartree self energy \f$ \Sigma_{ab}(\mathbf{k}) \f$ calculator

    Computes the Hartree self-energy of a static interaction as the sum

    \f[
        \Sigma_{ab}(\mathbf{k}) = \frac{1}{N_k}
          \sum_{\mathbf{q},cd} V_{abcd}(\mathbf{q}) \rho_{cd}(\mathbf{k} + \mathbf{q})
    \f]

    where \f$ \rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k}) \f$ is the density matrix of the
    single particle Green's function.

    @param v_k static interaction \f$ V_{abcd}(\mathbf{q}) \f$
    @param g_wk single particle Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$
    @return Hartree self-energy \f$ \Sigma_{ab}(\mathbf{k}) \f$
*/

  e_k_t hartree_sigma(chi_k_cvt v_k, g_wk_cvt g_wk);

  e_r_t hartree_sigma(chi_k_cvt v_k, e_r_cvt rho_r);
  
  /** Fock self energy \f$ \Sigma_{ab}(\mathbf{k}) \f$ calculator

    Computes the Fock self-energy of a static interaction as the sum

    \f[
        \Sigma_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{q}) \rho_{dc}(\mathbf{k} + \mathbf{q})
    \f]

    where \f$ \rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k}) \f$ is the density matrix of the
    single particle Green's function.

    @param v_k static interaction \f$ V_{abcd}(\mathbf{q}) \f$
    @param g_wk single particle Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$
    @return Fock self-energy \f$ \Sigma_{ab}(\mathbf{k}) \f$
*/

  e_k_t fock_sigma(chi_k_cvt v_k, g_wk_cvt g_wk);
  e_k_t fock_sigma(chi_k_cvt v_k, g_Dwk_cvt g_wk);

  e_r_t fock_sigma(chi_r_cvt v_r, e_r_cvt rho_r);

  /** Static GW self energy \f$ \Sigma_{ab}(\mathbf{k}) \f$ calculator

    Computes the static GW self-energy (equivalent to the Fock self-energy)
 
    @param v_k static interaction \f$ V_{abcd}(\mathbf{q}) \f$
    @param g_wk single particle Green's function \f$ G_{ab}(i\omega_n, \mathbf{k}) \f$
    @return Static GW self-energy (Fock) \f$ \Sigma_{ab}(\mathbf{k}) \f$
*/

  e_k_t gw_sigma(chi_k_cvt v_k, g_wk_cvt g_wk);
  e_k_t gw_sigma(chi_k_cvt v_k, g_Dwk_cvt g_wk);
  
  /** Dynamic GW self energy \f$ \Sigma(\tau, \mathbf{r}) \f$ calculator 

    Computes the GW self-energy as the product

    \f[
        \Sigma_{ab}(\tau, \mathbf{r}) =
          - \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})
    \f]

    @param W_tr interaction \f$ W_{abcd}(\tau, \mathbf{r}) \f$
    @param g_tr single particle Green's function \f$ G_{ab}(\tau, \mathbf{r}) \f$
    @return Dynamic GW self-energy \f$ \Sigma_{ab}(\tau, \mathbf{r}) \f$
 */

  g_tr_t gw_dynamic_sigma(chi_tr_cvt W_tr, g_tr_cvt g_tr);
  g_Dtr_t gw_dynamic_sigma(chi_Dtr_cvt W_tr, g_Dtr_cvt g_tr);

  /** some documentation */

  g_f_t g0w_dynamic_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, mesh::brzone::value_t kpoint);

  /** some documentation */

  g_fk_t g0w_dynamic_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, mesh::brzone kmesh);

  /** Real frequency GW self energy \f$ \Sigma(\omega, \mathbf{k}) \f$ calculator via the spectral representation

    Computes the spectral function of the dynamic part of the screened interaction
    
    \f[
        W^{(spec)}_{ab}(\omega, \mathbf{k}) = \frac{-1}{\pi} \text{Im}
          \left( W_{aabb}(\omega, \mathbf{k}) - V_{aabb}(\mathbf{k}) \right)
    \f]
          
    and constructs the dynamic part of the GW self energy via the spectral representation
    
    \f[
        \Sigma_{ab}(\omega, \mathbf{k}) = \frac{\delta_{\omega}}{N_k} \sum_{\mathbf{q}} \sum_{\omega'}
          U_{al}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{lb}(\mathbf{k}+\mathbf{q})
          W^{(spec)}_{ab}(\omega', \mathbf{q})
          \frac{n_B(\omega') + f(\epsilon_{\mathbf{k}+\mathbf{q}, l})}{\omega + i\delta + \omega' - \epsilon_{\mathbf{k}+\mathbf{q}, l} + \mu}
    \f]
          
    where \f$ \delta_{\omega} \f$ is the real-frequency mesh spacing and the \f$ U(\mathbf{k}) \f$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, i.e.

    \f[
       \sum_{\bar{a}b} U^\dagger_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}
    \f]
       
    @param mu chemical potential \f$ \mu \f$
    @param beta inverse temperature
    @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
    @param W_fk fully screened interaction \f$ W_{abcd}(\omega, \mathbf{k}) \f$
    @param v_k bare interaction \f$ V_{abcd}(\mathbf{k}) \f$
    @param delta broadening \f$ \delta \f$
    @return real frequency GW self-energy \f$ \Sigma_{ab}(\omega, \mathbf{k}) \f$
*/

  g_fk_t g0w_dynamic_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta);

  /** Some documentation */

  array<std::complex<double>, 2> g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k, mesh::brzone::value_t kpoint);

  /** Some documentation */

  e_k_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k, mesh::brzone kmesh);

  /** GW self energy \f$ \Sigma(\mathbf{k}) \f$ calculator for static interactions

    Computes the GW self-energy of a static interaction as the product

    \f[
        \Sigma_{ab}(\mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
          U_{al}(\mathbf{k}+\mathbf{q}) U^\dagger_{lb}(\mathbf{k}+\mathbf{q})
          V_{aabb}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l})
    \f]

    where the \f$ U(\mathbf{k}) \f$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, i.e.

    \f[
       \sum_{\bar{a}b} U^\dagger_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}
    \f]

    @param mu chemical potential \f$ \mu \f$
    @param beta inverse temperature
    @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
    @param v_k bare interaction \f$ V_{abcd}(\mathbf{k}) \f$
    @return static GW self-energy \f$ \Sigma_{ab}(\mathbf{k}) \f$
*/

  e_k_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k);

  /** some documentation */

  g_f_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, mesh::brzone::value_t kpoint);

  /** some documentation */

  g_fk_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta, mesh::brzone kmesh);

  /** Real frequency GW self energy \f$ \Sigma(\omega, \mathbf{k}) \f$ calculator via the spectral representation

    Computes the spectral function of the dynamic part of the screened interaction
    
    \f[
        W^{(spec)}_{ab}(\omega, \mathbf{k}) = \frac{-1}{\pi} \text{Im}
          \left( W_{aabb}(\omega, \mathbf{k}) - V_{aabb}(\mathbf{k}) \right)
    \f]
          
    and constructs the GW self energy via the spectral representation
    
    \f[
        \Sigma_{ab}(\omega, \mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
          U_{al}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{lb}(\mathbf{k}+\mathbf{q})
          V_{aabb}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l}) \\
        + \frac{\delta_{\omega}}{N_k} \sum_{\mathbf{q}} \sum_{\omega'}
          U_{al}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{lb}(\mathbf{k}+\mathbf{q})
          W^{(spec)}_{ab}(\omega', \mathbf{q})
          \frac{n_B(\omega') + f(\epsilon_{\mathbf{k}+\mathbf{q}, l})}{\omega + i\delta + \omega' - \epsilon_{\mathbf{k}+\mathbf{q}, l} + \mu}
    \f]
          
    where \f$ \delta_{\omega} \f$ is the real-frequency mesh spacing and the \f$ U(\mathbf{k}) \f$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$, i.e.

    \f[
       \sum_{\bar{a}b} U^\dagger_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}
    \f]
       
    @param mu chemical potential \f$ \mu \f$
    @param beta inverse temperature
    @param e_k discretized lattice dispersion \f$ \epsilon_{\bar{a}b}(\mathbf{k}) \f$
    @param W_fk fully screened interaction \f$ W_{abcd}(\omega, \mathbf{k}) \f$
    @param v_k bare interaction \f$ V_{abcd}(\mathbf{k}) \f$
    @param delta broadening \f$ \delta \f$
    @return real frequency GW self-energy \f$ \Sigma_{ab}(\omega, \mathbf{k}) \f$
*/

  g_fk_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta);

} // namespace triqs_tprf
