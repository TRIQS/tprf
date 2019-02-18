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
	 \Gamma_{A\bar{a}\bar{b}B}(\mathbf{k}-\mathbf{k}', i\nu - i\nu')
	 \\ \times
	 G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')

     @param chi_pp particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}\bar{c}d}(\mathbf{k}, i\nu_n)`
     @param g_kw single particle Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\nu_n)`
     @param delta_kw pairing self-energy :math:`\Delta_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n)`
     @return Gives the result of the product :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG \Delta`

  */

  gk_iw_t eliashberg_product(chi_wk_vt Gamma_pp, gk_iw_vt g_wk, gk_iw_vt delta_wk);

}
