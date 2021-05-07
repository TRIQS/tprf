/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2017 by M. Ferrero, O. Parcollet
 * Copyright (C) 2018- by Simons Foundation
 *               authors : O. Parcollet, N. Wentzell, H. U.R. Strand
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

namespace triqs_tprf::fourier {

  using fourier_plan = std::unique_ptr<void, void (*)(void *)>;

  fourier_plan _fourier_base_plan(nda::array_const_view<dcomplex, 2> in, nda::array_const_view<dcomplex, 2> out, int rank, int *dims, int fftw_count,
                                  int fftw_backward_forward);

  void _fourier_base(nda::array_const_view<dcomplex, 2> in, nda::array_view<dcomplex, 2> out, fourier_plan &p);

} // namespace triqs_tprf::fourier
