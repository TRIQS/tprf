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

/** Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator 
    for static momentum-dependent interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) = 
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

    Instead of returning :math:`W^{(full)}` we return the dynamical/retarded part
    :math:`W^{(r)}` (with zero high-frequency offset)
    
    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) = 
            W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) - V_{abcd}(\mathbf{k})

    @param PI_wk polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_k static interaction :math:`V_{abcd}(\mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */

chi_wk_t dynamical_screened_interaction_W_wk(chi_wk_cvt PI_wk, chi_k_cvt V_k);

/** Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator 
    for static momentum-dependent interactions :math:`V(\mathbf{k})` and 
    known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) = 
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          V_{hgcd}(\mathbf{k})

    Instead of returning :math:`W^{(full)}` we return the dynamical/retarded part
    :math:`W^{(r)}` (with zero high-frequency offset)
    
    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) = 
            W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) - V_{abcd}(\mathbf{k})

    @param chi_wk polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_k static interaction :math:`V_{abcd}(\mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */

chi_wk_t dynamical_screened_interaction_W_wk_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k);

/** GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator 

    Fourier transforms the screened interaction and the single-particle
    Green's function to imagiary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W^{(r)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{(r)}_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W^{(r)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma_{ab}(\tau, \mathbf{r}) \right\}

    @param V_k static bare interaction :math:`V_{abcd}(\mathbf{k})`
    @param Wr_wk retarded screened interaction :math:`W^{(r)}_{abcd}(i\omega_n, \mathbf{k})`
    @param g_wk single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`
    @return GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`
 */

g_wk_t gw_sigma_wk_serial_fft(chi_wk_cvt Wr_wk, g_wk_cvt g_wk);

/** GW self energy :math:`\Sigma(\tau, \mathbf{r})` calculator 

    Computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W^{(r)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    @param Wr_tr retarded screened interaction :math:`W^{(r)}_{abcd}(\tau, \mathbf{r})`
    @param g_tr single particle Green's function :math:`G_{ab}(\tau, \mathbf{r})`
    @return GW self-energy :math:`\Sigma_{ab}(\tau, \mathbf{r})`
 */

g_tr_t gw_sigma_tr(chi_tr_cvt Wr_tr, g_tr_cvt g_tr);

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
    @param fmesh real frequency mesh
    @param Wr_fk fully screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
    @param V_k bare interaction :math:`V_{abcd}(\mathbf{k})`
    @param delta broadening :math:`\delta`
    @return real frequency GW self-energy :math:`\Sigma_{ab}(\omega, \mathbf{k})`
*/

g_fk_t gw_sigma_fk_g0w0_spectral(double mu, double beta, e_k_cvt e_k, gf_mesh<refreq> mesh, chi_fk_cvt Wr_fk, chi_k_cvt v_k, double delta);


} // namespace triqs_tprf
