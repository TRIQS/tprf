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

namespace tprf {

/** Screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator 
    for static momentum dependent interaction :math:`V(\mathbf{k})`

    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) = 
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          W_{hgcd}(i\omega_n, \mathbf{k})

    @param PI_wk polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_k static interaction :math:`V_{abcd}(\mathbf{k})`
    @return screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */

chi_wk_t screened_interaction_W(chi_wk_vt PI_wk, chi_k_vt V_k);

/** GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator 

    Fourier transforms the screened interaction and the single-particle
    Green's function to imagiary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma_{ab}(\tau, \mathbf{r}) \right\}

    @param W_wk screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
    @param g_wk single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`
    @return GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`
 */

g_wk_t gw_self_energy(chi_wk_vt W_wk, g_wk_vt g_wk);

} // namespace tprf
