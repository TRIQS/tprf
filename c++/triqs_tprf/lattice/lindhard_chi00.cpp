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

#include <triqs/arrays/linalg/eigenelements.hpp>

#include "common.hpp"
#include "lindhard_chi00.hpp"

namespace triqs_tprf {

// ----------------------------------------------------
// chi00 bubble in analytic form

double fermi(double e) { return 1. / (exp(e) + 1.); }

chi_wk_t lindhard_chi00_wk(e_k_cvt e_k, int nw,
                          double beta, double mu) {

  auto kmesh = e_k.mesh();
  int nb = e_k.target().shape()[0];

  chi_wk_t chi_wk{{{beta, Boson, nw}, kmesh}, {nb, nb, nb, nb}};
  for (auto const & [ w, k ] : chi_wk.mesh())
    chi_wk[w, k] = 0.;

  auto wmesh = std::get<0>(chi_wk.mesh());

  for (auto const &k : kmesh) {

    //std::cout << "kidx, k = " << k.linear_index() << ", " << k << "\n";
	
    //for (auto const &q : kmesh) { // can not do range-based for loops with OpenMP

#pragma omp parallel for 
    for (unsigned int qidx = 0; qidx < kmesh.size(); qidx++) {
      auto q_iter = kmesh.begin();
      q_iter += qidx;
      auto q = *q_iter;

      // -- If this is moved out to the k-loop the threading breaks?!?
      matrix<std::complex<double>> e_k_mat(e_k[k] - mu);
      auto eig_k = linalg::eigenelements(e_k_mat);
      auto ek = eig_k.first;
      auto Uk = eig_k.second;

      matrix<std::complex<double>> e_kq_mat(e_k(k + q) - mu);
      auto eig_kq = linalg::eigenelements(e_kq_mat);
      auto ekq = eig_kq.first;
      auto Ukq = eig_kq.second;

      for (int i : range(nb)) {
        for (int j : range(nb)) {

          double de = ekq(j) - ek(i);
          double dn = fermi(ek(i) * beta) - fermi(ekq(j) * beta);

          for (auto const &w : wmesh) {

            std::complex<double> total_factor = dn / (w + de);

            double tol = 1e-10;

            if (abs(std::complex<double>(w)) < tol &&
                abs(de) < tol) { // w=0, de=0, 2nd order pole

              // -- analytic first derivative of the fermi distribution function
              // -- evaluated at ek(i)

              double cosh_be = cosh(0.5 * beta * ek(i));
              total_factor = beta / (4. * cosh_be * cosh_be);
            }

            chi_wk[w, q](a, b, c, d)
                << chi_wk[w, q](a, b, c, d) + Uk(i, a) * dagger(Uk)(d, i) *
                                               Ukq(j, c) * dagger(Ukq)(b, j) *
                                               total_factor;
          } // w
        }   // j
      }     // i
    }       // q
  }         // k

  chi_wk /= kmesh.size();

  return chi_wk;
} 

} // namespace triqs_tprf
