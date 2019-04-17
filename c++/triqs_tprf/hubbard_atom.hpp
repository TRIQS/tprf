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

namespace triqs_tprf {
namespace hubbard_atom {

/** Single-particle Green's function of the Hubbard atom at half-filling

     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
     please cite the paper if you use this function!

     .. math::
         G(i\omega_n) = \frac{1}{i\omega_n - \frac{U^2}{4 i\omega_n}}

     @param nw number of Matsubara frequencies
     @param beta inverse temperature
     @param U Hubbard U interaction parmeter
     @return g single-particle Green's function of the Hubbard atom :math:`G(i\omega_n)`
  */
g_iw_t single_particle_greens_function(int nw, double beta, double U);

/** Magnetic susceptibility of the Hubbard atom at half-filling :math:`\chi(\omega, \nu, \nu')`

     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
     please cite the paper if you use this function!

     @param nw number of bosonic Matsubara frequencies
     @param nwf number of fermionic Matsubara frequencies
     @param beta inverse temperature
     @param U Hubbard U interaction parmeter
     @return chi magnetic susceptibility :math:`\chi(\omega, \nu, \nu')`
  */
g2_iw_t chi_ph_magnetic(int nw, int nwf, double beta, double U);

/** Magnetic vertex function in the particle-hole channel of the Hubbard atom at half-filling :math:`\Gamma(\omega, \nu, \nu')`

    Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
    please cite the paper if you use this function!

    @param nw number of bosonic Matsubara frequencies
    @param nwf number of fermionic Matsubara frequencies
    @param beta inverse temperature
    @param U Hubbard U interaction parmeter
    @return gamma magnetic susceptibility :math:`\Gamma(\omega, \nu, \nu')`
  */
g2_iw_t gamma_ph_magnetic(int nw, int nwf, double beta, double U);

} // namespace hubbard_atom
} // namespace triqs_tprf
