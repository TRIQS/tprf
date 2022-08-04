/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, The Simons Foundation
 * Authors: H. U.R. Strand, Y. in 't Veld, N. Wentzell
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

  /** Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) = 
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

    @param PI_wk polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_k static bare interaction :math:`V_{abcd}(\mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */
  chi_wk_t dynamical_screened_interaction_W(chi_wk_cvt PI_wk, chi_k_cvt V_k);

  /** Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) = 
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \Pi_{fegh}(\omega, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(\omega, \mathbf{k})

    @param PI_fk polarization bubble :math:`\Pi_{abcd}(\omega, \mathbf{k})`
    @param V_k static bare interaction :math:`V_{abcd}(\mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
 */
  chi_fk_t dynamical_screened_interaction_W(chi_fk_cvt PI_fk, chi_k_cvt V_k);

  /** Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(i\omega_n, \mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) = 
        V_{abcd}(i\omega_n, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(i\omega_n, \mathbf{k}) \cdot
        \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

    @param PI_wk polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_wk dynamic bare interaction :math:`V_{abcd}(i\omega_n, \mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */
  chi_wk_t dynamical_screened_interaction_W(chi_wk_cvt PI_wk, chi_wk_cvt V_wk);

  /** Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(\omega, \mathbf{k})`.

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) = 
        V_{abcd}(\omega, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(\omega, \mathbf{k}) \cdot
        \Pi_{fegh}(\omega, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(\omega, \mathbf{k})

    @param PI_fk polarization bubble :math:`\Pi_{abcd}(\omega, \mathbf{k})`
    @param V_fk dynamic bare interaction :math:`V_{abcd}(\omega, \mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
 */
  chi_fk_t dynamical_screened_interaction_W(chi_fk_cvt PI_fk, chi_fk_cvt V_fk);

  /** Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})` and known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) = 
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        V_{hgcd}(\mathbf{k})

    @param chi_wk generalized susceptibility :math:`\chi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_k static bare interaction :math:`V_{abcd}(\mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */
  chi_wk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_k_cvt V_k);

  /** Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})` and known generalized susceptibility :math:`\chi(\omega, \mathbf{k})`

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) = 
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \chi_{fegh}(\omega, \mathbf{k}) \cdot
        V_{hgcd}(\mathbf{k})

    @param chi_fk generalized susceptibility :math:`\chi_{abcd}(\omega, \mathbf{k})`
    @param V_k static bare interaction :math:`V_{abcd}(\mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
 */
  chi_fk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk_cvt chi_fk, chi_k_cvt V_k);

  /** Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(i\omega_n, \mathbf{k})` and known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) = 
        V_{abcd}(i\omega_n, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(i\omega_n, \mathbf{k}) \cdot
        \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        V_{hgcd}(i\omega_n, \mathbf{k})

    @param chi_wk generalized susceptibility :math:`\chi_{abcd}(i\omega_n, \mathbf{k})`
    @param V_wk dynamic bare interaction :math:`V_{abcd}(i\omega_n, \mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
 */
  chi_wk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk_cvt chi_wk, chi_wk_cvt V_wk);

  /** Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(\omega, \mathbf{k})` and known generalized susceptibility :math:`\chi(\omega, \mathbf{k})`

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) = 
        V_{abcd}(\omega, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(\omega, \mathbf{k}) \cdot
        \chi_{fegh}(\omega, \mathbf{k}) \cdot
        V_{hgcd}(\omega, \mathbf{k})

    @param chi_fk generalized susceptibility :math:`\chi_{abcd}(\omega, \mathbf{k})`
    @param V_fk dynamic bare interaction :math:`V_{abcd}(\omega, \mathbf{k})`
    @return dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
 */
  chi_fk_t dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk_cvt chi_fk, chi_fk_cvt V_fk);

} // namespace triqs_tprf
