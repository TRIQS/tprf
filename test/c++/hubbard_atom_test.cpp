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

#include <triqs/test_tools/gfs.hpp>

#include "hubbard_atom.hpp"

using namespace tprf;
using namespace tprf::hubbard_atom;

TEST(hubbard_atom, single_particle_greens_function) {
  double eps = 1e-9;
  double beta = 2.0;
  double U = 5.0;
  int nw = 10000;
  auto G_iw = single_particle_greens_function(nw, beta, U);

  for( auto const &iw : G_iw.mesh() ) {
    auto ref_val = 1./(iw - U*U/(4.*iw));
    EXPECT_NEAR(G_iw[iw](0,0).real(), ref_val.real(), eps);
    EXPECT_NEAR(G_iw[iw](0,0).imag(), ref_val.imag(), eps);
  }
}

TEST(hubbard_atom, chi_ph_magnetic) {
  double eps = 1e-9;
  double beta = 2.234;
  double U = 5.0;
  int nw = 4;
  int nwf = 10;
  auto chi = chi_ph_magnetic(nw, nwf, beta, U);

  std::cerr << "Lacking reference data FIXME!\n";
}

MAKE_MAIN;
