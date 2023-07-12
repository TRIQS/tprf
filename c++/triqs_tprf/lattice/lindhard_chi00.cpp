/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, H. U.R. Strand, Y. in 't Veld, M. Rösner
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

#include <nda/linalg.hpp>

#include "common.hpp"
#include "lindhard_chi00.hpp"
#include "lattice_utility.hpp"
#include "../mpi.hpp"

namespace triqs_tprf {

  // ----------------------------------------------------
  // chi00 bubble in analytic form

  template<typename chi_t, typename mesh_t>
  chi_t lindhard_chi00_template(e_k_cvt e_k, mesh_t mesh, double beta, double mu, double delta=0.) {

    auto wmesh = mesh;
    auto kmesh = e_k.mesh();
    int nb     = e_k.target().shape()[0];
    std::complex<double> idelta(0.0, delta);

    chi_t chi_wk{{wmesh, kmesh}, {nb, nb, nb, nb}};
    for (auto [w, k] : chi_wk.mesh()) chi_wk[w, k] = 0.;

    auto arr = mpi_view(kmesh);

#pragma omp parallel for
    for (unsigned int qidx = 0; qidx < kmesh.size(); qidx++) {
      auto q = *std::next(kmesh.begin(), qidx);

      for (auto k : arr) {

        // -- If this is moved out to the k-loop the threading breaks?!?
        matrix<std::complex<double>> e_k_mat(e_k[k] - mu);
        auto [ek, Uk] = linalg::eigenelements(e_k_mat);

        matrix<std::complex<double>> e_kq_mat(e_k(k + q) - mu);
        auto [ekq, Ukq] = linalg::eigenelements(e_kq_mat);

        for (int i : range(nb)) {
          for (int j : range(nb)) {

            double de = ekq(j) - ek(i);
            double dn = fermi(ek(i) * beta) - fermi(ekq(j) * beta);

            for (auto w : wmesh) {

              std::complex<double> total_factor;

              double tol = 1e-10;
              if (abs(std::complex<double>(w) + idelta) < tol && abs(de) < tol) {
                // w=0, de=0, 2nd order pole

                // -- analytic first derivative of the fermi distribution function
                // -- evaluated at ek(i)

                double cosh_be = cosh(0.5 * beta * ek(i));
                total_factor   = beta / (4. * cosh_be * cosh_be);
              } else {
                total_factor = dn / (w + idelta + de);
              }

              chi_wk[w, q](a, b, c, d) << chi_wk[w, q](a, b, c, d) + Uk(a, i) * dagger(Uk)(i, d) * Ukq(c, j) * dagger(Ukq)(j, b) * total_factor;
            } // w
          }   // j
        }     // i
      }       // q
    }         // k

    chi_wk = mpi::all_reduce(chi_wk);
    chi_wk /= kmesh.size();

    return chi_wk;
  }

  chi_wk_t lindhard_chi00(e_k_cvt e_k, mesh::imfreq mesh, double mu) {
    if (mesh.statistic() != Boson) TRIQS_RUNTIME_ERROR << "lindhard_chi00: statistic is incorrect.\n";
    return lindhard_chi00_template<chi_wk_t, mesh::imfreq>(e_k, mesh, mesh.beta(), mu);
  }

  chi_Dwk_t lindhard_chi00(e_k_cvt e_k, mesh::dlr_imfreq mesh, double mu) {
    if (mesh.statistic() != Boson) TRIQS_RUNTIME_ERROR << "lindhard_chi00: statistic is incorrect.\n";
    return lindhard_chi00_template<chi_Dwk_t, mesh::dlr_imfreq>(e_k, mesh, mesh.beta(), mu);
  }

  chi_fk_t lindhard_chi00(e_k_cvt e_k, mesh::refreq mesh, double beta, double mu, double delta) {
    return lindhard_chi00_template<chi_fk_t, mesh::refreq>(e_k, mesh, beta, mu, delta);
  }
  
} // namespace triqs_tprf
