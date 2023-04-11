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

#include "fourier_interpolation.hpp"

namespace triqs_tprf {

// ----------------------------------------------------
// fourier interpolation

array<std::complex<double>, 6> cluster_mesh_fourier_interpolation(array<double, 2> k_vecs, chi_wr_cvt chi) {

  int nk = k_vecs.shape()[0];
  int nb = chi.target().shape()[0];
  int nw = std::get<0>(chi.mesh()).size();

  auto wmesh = std::get<0>(chi.mesh());
  auto rmesh = std::get<1>(chi.mesh());

  array<std::complex<double>, 6> chi_out(nw, nk, nb, nb, nb, nb);

#pragma omp parallel for
  for(int kidx = 0; kidx < nk; kidx++) {
    
    chi_out(range::all, kidx, range::all, range::all, range::all, range::all) *= 0.;

    auto k = k_vecs(kidx, range::all);
    
    for (auto const &r : rmesh ) {

    /*
#pragma omp parallel for
    for (int idx = 0; idx < rmesh.size(); idx++) {
      auto iter = rmesh.begin(); iter += idx; auto r = *iter;
    */

      auto dot_prod = k(0)*r(0) + k(1)*r(1) + k(2)*r(2);
      auto exponent = exp( - std::complex<double>(0., dot_prod) );

      for (auto const &w : wmesh) {
	int widx = w.linear_index();
	for( int a : range(nb) )
	for( int b : range(nb) )
	for( int c : range(nb) )
	for( int d : range(nb) )
	  chi_out(widx, kidx, a, b, c, d) += exponent * chi[w, r](a, b, c, d);
      }
    }
  }
  return chi_out;
}

} // namespace triqs_tprf
