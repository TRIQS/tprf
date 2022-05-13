/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2019, The Simons Foundation
 * Authors: H. U.R. Strand
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

/** GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator for dynamic interactions

    Splits the interaction into a dynamic and a static part
    
    .. math ::
        W_{abcd}(i\omega_n, \mathbf{k}) = 
            W^{dyn}_{abcd}(i\omega_n, \mathbf{k})
            + V_{abcd}(\mathbf{k})

    by fitting the high-frequency tail.

    Fourier transforms the dynamic part of the interaction and the 
    single-particle Green's function to imaginary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{dyn}_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma^{dyn}_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma_{ab}(\tau, \mathbf{r}) \right\}

    The self-energy of the static part of the interaction is calculated
    as the sum

    .. math::
        \Sigma^{stat}_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{mathbf{q}} V_{abab}(\mathbf{k}) \rho(G(i\omega_n, \mathbf{k+q}))_{ab}

    where :math:`\rho(G(i\omega_n, \mathbf{k+q}))` is the density matrix of the
    single particle Green's function.

    The total GW self-energy is given by

    ,, math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) = 
          \Sigma^{dyn}_{ab}(i\omega_n, \mathbf{k})
          + \Sigma^{stat}_{ab}(\mathbf{k})

    @param W_wk interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
    @param g_wk single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`
    @return GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`
 */

g_wk_t gw_sigma(chi_wk_cvt W_wk, g_wk_cvt g_wk);

/** GW self energy :math:`\Sigma(\mathbf{k})` calculator for static interactions

    Computes the GW self-energy of a static interaction as the sum

    .. math::
        \Sigma_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{mathbf{q}} V_{abab}(\mathbf{k}) \rho(G(i\omega_n, \mathbf{k+q}))_{ab}

    where :math:`\rho(G(i\omega_n, \mathbf{k+q}))` is the density matrix of the
    single particle Green's function.

    @param V_k static interaction :math:`V_{abcd}(\mathbf{k})`
    @param g_wk single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`
    @return GW self-energy :math:`\Sigma_{ab}(\mathbf{k})`
*/

e_k_t gw_sigma(chi_k_cvt v_k, g_wk_cvt g_wk);

/** GW self energy :math:`\Sigma(\tau, \mathbf{r})` calculator 

    Computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    @param W_tr interaction :math:`W_{abcd}(\tau, \mathbf{r})`
    @param g_tr single particle Green's function :math:`G_{ab}(\tau, \mathbf{r})`
    @return GW self-energy :math:`\Sigma_{ab}(\tau, \mathbf{r})`
 */

g_tr_t gw_sigma(chi_tr_cvt W_tr, g_tr_cvt g_tr);


/** Real frequency GW self energy :math:`\Sigma(\omega, \mathbf{k})` calculator via the spectral representation

    Computes the spectral function of the dynamic part of the screened interaction
    
    .. math::
        W^{(spec)}_{ab}(\omega, \mathbf{k}) = \frac{-1}{\pi} \text{Im}
          \left( W_{abab}(\omega, \mathbf{k}) - V_{abab}(\mathbf{k}) \right)
          
    and constructs the GW self energy via the spectral representation
    
    .. math::
        \Sigma_{ab}(\omega, \mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
          U_{la}(\mathbf{k}+\mathbf{q}) U^\dagger_{bl}(\mathbf{k}+\mathbf{q})
          V_{abab}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l})
        + \frac{\delta_{\omega}}{N_k} \sum_{\mathbf{q}} \sum_{\omega^'}
          U_{la}(\mathbf{k}+\mathbf{q}) U^\dagger_{bl}(\mathbf{k}+\mathbf{q})
          W^{(spec)}_{ab}(\omega^', \mathbf{q})
          \frac{n_B(\omega^') + f(\epsilon_{\mathbf{k}+\mathbf{q}, l})}
          {\omega + i\delta + \omega^' - \epsilon_{\mathbf{k}+\mathbf{q}, l}) + \mu}
          
    where the $U(\mathbf{k})$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation $\epsilon_{\bar{a}b}(\mathbf{k})$, i.e.

    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}
       
    @param mu chemical potential :math:`\mu`
    @param beta inverse temperature
    @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
    @param mesh real frequency mesh
    @param W_fk fully screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
    @param V_k bare interaction :math:`V_{abcd}(\mathbf{k})`
    @param delta broadening :math:`\delta`
    @return real frequency GW self-energy :math:`\Sigma_{ab}(\omega, \mathbf{k})`
*/

g_fk_t g0w_sigma(double mu, double beta, e_k_cvt e_k, gf_mesh<refreq> mesh, chi_fk_cvt W_fk, chi_k_cvt v_k, double delta);

/** GW self energy :math:`\Sigma(\mathbf{k})` calculator for static interactions

    Computes the GW self-energy of a static interaction as the product

    .. math::
        \sigma_{ab}(\mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
          U_{la}(\mathbf{k}+\mathbf{q}) U^\dagger_{bl}(\mathbf{k}+\mathbf{q})
          V_{abab}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l})

    where the $U(\mathbf{k})$ matrices are the diagonalizing unitary transform of the matrix valued 
    dispersion relation $\epsilon_{\bar{a}b}(\mathbf{k})$, i.e.

    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}

    @param mu chemical potential :math:`\mu`
    @param beta inverse temperature
    @param e_k discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
    @param V_k bare interaction :math:`V_{abcd}(\mathbf{k})`
    @return static GW self-energy :math:`\Sigma_{ab}(\mathbf{k})`
*/

e_k_t g0w_sigma(double mu, double beta, e_k_cvt e_k, chi_k_cvt v_k);

} // namespace triqs_tprf
