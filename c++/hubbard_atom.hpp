/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
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

#include "types.hpp"

namespace tprf {
namespace hubbard_atom {
  
  typedef std::complex<double> val_t; 
  typedef gf<imfreq, tensor_valued<4>> temp_1d_t;
  typedef gf<cartesian_product<imfreq, imfreq>, tensor_valued<4>> temp_2d_t;

  g_iw_t single_particle_greens_function(int nw, double beta, double U);
  g2_iw_t chi_ph_magnetic(int nw, int nwf, double beta, double U);
  
} // namespace hubbard_atom
} // namespace tprf
