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
  /** Split mtrix-valued Green's function into dynamic and constant parts by tail fitting

  @param chi_wk : general matrix-valued Green's function :math:`\chi(i\omega_n, \mathbf{k})`. 
  @return Tuple of chi_dyn_wk, the dynamic part of chi, which converges to zero for :math:`\omega_n \rightarrow \infty`, and chi_const_k, the part of chi that is constant in Matsubara frequency space :math:`\chi(\mathbf{k})`.
  */
  std::tuple<chi_wk_t, chi_k_t> split_into_dynamic_wk_and_constant_k(chi_wk_cvt chi_wk);
}
