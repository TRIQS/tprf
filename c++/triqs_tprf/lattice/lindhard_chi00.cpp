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
#include "../mpi.hpp"

namespace triqs_tprf {

// ----------------------------------------------------
// helper functions

double fermi(double e) { return 1. / (exp(e) + 1.); }

// ----------------------------------------------------
// chi00 bubble in analytic form in real frequencies

chi_fk_t lindhard_chi00_fk(h_k_cvt h_k, gf_mesh<refreq> mesh, double beta, 
                           double mu, double delta) {

  auto kmesh = h_k.mesh();
  int nb = h_k.target().shape()[0];

  chi_fk_t chi_fk{{mesh, kmesh}, {nb, nb, nb, nb}};
  
  auto fmesh = std::get<0>(chi_fk.mesh());
  
  for (auto const & [ f, k ] : chi_fk.mesh())
    chi_fk[f, k] = 0.;
    
  std::complex<double> idelta(0.0, delta);
    
  auto arr = mpi_view(kmesh);

  #pragma omp parallel for
  for (int qidx = 0; qidx < kmesh.size(); qidx++) {
  
      auto q_iter = kmesh.begin();
      q_iter += qidx;
      auto q = *q_iter;
  
      for (int kidx = 0; kidx < arr.size(); kidx++) {
      
          auto k = arr(kidx);
          
          matrix<std::complex<double>> h_k_mat(h_k(k));
          auto eig_k = linalg::eigenelements(h_k_mat);
          auto ek = eig_k.first;
          auto Uk = eig_k.second;
          
          matrix<std::complex<double>> h_kq_mat(h_k(k + q));
          auto eig_kq = linalg::eigenelements(h_kq_mat);
          auto ekq = eig_kq.first;
          auto Ukq = eig_kq.second;
          
          for (int i : range(nb)) {
            for (int j : range(nb)) {
          
              double dn = fermi((ek(i) - mu) * beta) - fermi((ekq(j) - mu) * beta);
              double de = ekq(j) - ek(i);
              
              for (auto const &f : fmesh) {
      
                chi_fk[f, q](a, b, c, d) 
                    << chi_fk[f, q](a, b, c, d) + Uk( i, a) * dagger(Uk )(d, i) 
                                                * Ukq(j, c) * dagger(Ukq)(b, j)
                                                * dn / (f + idelta + de);
                  
              } // f
              
          } // i
        } // j
          
      } // k (mpi)
  
  } // q (omp)
  
  chi_fk = mpi::all_reduce(chi_fk);
  
  chi_fk /= kmesh.size();
    
  return chi_fk;

} 

// ----------------------------------------------------
// chi00 bubble in analytic form

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
    for (int qidx = 0; qidx < kmesh.size(); qidx++) {
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
