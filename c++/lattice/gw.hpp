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
namespace gw {

/** W calculator

    Computes

    .. math::
        W(i\omega_n) = ...

    @param X math parameter :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{k},
   i\nu_n)`
    @return Gives the result of :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG
   \Delta`

 */

chi_wk_t screened_interaction_W(chi_wk_vt PI_wk, chi_k_vt V_k);
  
} // namespace gw
} // namespace tprf
