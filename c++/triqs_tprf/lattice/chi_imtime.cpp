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

#include "common.hpp"
#include "../mpi.hpp"
#include "chi_imtime.hpp"

#include "../fourier/fourier.hpp"
#include "fourier.hpp"

namespace triqs_tprf {

  namespace {
    using namespace fourier;
  }

// ----------------------------------------------------
// chi0 bubble in imaginary time

chi_tr_t chi0_tr_from_grt_PH(g_tr_cvt g_tr) {

  auto _ = all_t{};

  auto tmesh = std::get<0>(g_tr.mesh());
  auto rmesh = std::get<1>(g_tr.mesh());

  int nb = g_tr.target().shape()[0];
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_tr_t chi0_tr{{{beta, Boson, ntau}, rmesh}, {nb, nb, nb, nb}};

  auto g_target = g_tr.target();
  auto chi_target = chi0_tr.target();
  
  // -- This does not work on the boundaries!! The eval wraps to the other
  // regime!
  // -- gt(beta) == gt(beta + 0^+)
  // chi0_tr(tau, r)(a, b, c, d) << g_tr(tau, r)(d, a) * g_tr(-tau, -r)(b, c);

  //for (auto const &r : rmesh) {

  auto arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto & r = arr(idx);

    auto chi0_t = make_gf<imtime>({beta, Boson, ntau}, chi_target);
    auto g_pr_t = make_gf<imtime>(tmesh, g_target);
    auto g_mr_t = make_gf<imtime>(tmesh, g_target);

#pragma omp critical
    {
      g_pr_t = g_tr[_, r];
      g_mr_t = g_tr[_, -r];
    }
    
    for (auto const &t : tmesh)
      chi0_t[t](a, b, c, d) << g_pr_t(t)(d, a) * g_mr_t(beta - t)(b, c);

#pragma omp critical
    chi0_tr[_, r] = chi0_t;
    
  }

  chi0_tr = mpi_all_reduce(chi0_tr);
  return chi0_tr;
}

// -- optimized version for w=0
chi_wr_t chi0_w0r_from_grt_PH(g_tr_cvt g_tr) {

  auto _ = all_t{};

  auto tmesh = std::get<0>(g_tr.mesh());
  auto rmesh = std::get<1>(g_tr.mesh());

  int nw = 1;
  int nb = g_tr.target().shape()[0];
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_wr_t chi0_wr{{{beta, Boson, nw}, rmesh}, {nb, nb, nb, nb}};

  auto g_target = g_tr.target();
  auto chi_target = chi0_wr.target();
  
  auto arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto & r = arr(idx);

    auto chi0_t = make_gf<imtime>({beta, Boson, ntau}, chi_target);
    auto g_pr_t = make_gf<imtime>(tmesh, g_target);
    auto g_mr_t = make_gf<imtime>(tmesh, g_target);

#pragma omp critical
    {
      g_pr_t = g_tr[_, r];
      g_mr_t = g_tr[_, -r];
    }
    
    for (auto const &t : tmesh)
      chi0_t[t](a, b, c, d) << g_pr_t(t)(d, a) * g_mr_t(beta - t)(b, c);

    auto int_chi0 = chi_trapz_tau(chi0_t);
    
#pragma omp critical
    chi0_wr[0, r] = int_chi0;

  }

  chi0_wr = mpi_all_reduce(chi0_wr);
  return chi0_wr;
}  

chi_t_t::zero_t chi_trapz_tau(chi_t_cvt chi_t) {

  auto tmesh = chi_t.mesh();
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_t_t::zero_t I = chi_t.get_zero();
    
  // -- Trapetzoidal integration

  for (auto const &t : tmesh) I += chi_t[t];

  auto boundary_iter = tmesh.begin();
  I -= 0.5 * chi_t[*boundary_iter];

  boundary_iter += ntau - 1;
  I -= 0.5 * chi_t[*boundary_iter];

  I *= beta / (ntau-1);
  return I;  
}

// -- specialized calc for w=0
chi_wr_t chi_w0r_from_chi_tr(chi_tr_cvt chi_tr) {

  int nb = chi_tr.target().shape()[0];

  auto tmesh = std::get<0>(chi_tr.mesh());
  auto rmesh = std::get<1>(chi_tr.mesh());

  int nw = 1;
  double beta = tmesh.domain().beta;

  chi_wr_t chi_wr{{{beta, Boson, nw}, rmesh}, {nb, nb, nb, nb}};

  auto arr = mpi_view(rmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto & r = arr(idx);

    auto _ = all_t{};
    auto I = chi_trapz_tau(chi_tr[_, r]);

#pragma omp critical
    chi_wr[0, r] = I;
  }

  chi_wr = mpi_all_reduce(chi_wr);
  return chi_wr;
}

chi_wr_t chi_wr_from_chi_tr(chi_tr_cvt chi_tr, int nw) {
  auto chi_wr = fourier_tr_to_wr_general_target(chi_tr, nw);
  return chi_wr;
}

chi_tr_t chi_tr_from_chi_wr(chi_wr_cvt chi_wr, int ntau) {
  auto chi_tr = fourier_wr_to_tr_general_target(chi_wr, ntau);
  return chi_tr;
}  

chi_wk_t chi_wk_from_chi_wr(chi_wr_cvt chi_wr) {
  auto chi_wk = fourier_wr_to_wk_general_target(chi_wr);
  return chi_wk;
}

chi_wr_t chi_wr_from_chi_wk(chi_wk_cvt chi_wk) {
  auto chi_wr = fourier_wk_to_wr_general_target(chi_wk);  
  return chi_wr;
}  
  
  /*
chi_wk_t chi_wk_from_chi_wr(chi_wr_cvt chi_wr) {

  // auto target = chi_wr.target();
  int nb = chi_wr.target().shape()[0];

  auto wmesh = std::get<0>(chi_wr.mesh());
  auto rmesh = std::get<1>(chi_wr.mesh());

  auto kmesh = gf_mesh<brillouin_zone>{brillouin_zone{rmesh.domain()}, rmesh.periodization_matrix};
  
  chi_wk_t chi_wk{{wmesh, kmesh}, {nb, nb, nb, nb}};

  auto _ = all_t{};
  for (auto const &w : wmesh)
    chi_wk[w, _] = triqs::gfs::fourier(chi_wr[w, _]);

  return chi_wk;
}

chi_wr_t chi_wr_from_chi_wk(chi_wk_cvt chi_wk) {

  int nb = chi_wk.target().shape()[0];

  auto wmesh = std::get<0>(chi_wk.mesh());
  auto kmesh = std::get<1>(chi_wk.mesh());

  auto rmesh = gf_mesh<cyclic_lattice>{bravais_lattice{kmesh.domain()}, kmesh.periodization_matrix};
  
  chi_wr_t chi_wr{{wmesh, rmesh}, {nb, nb, nb, nb}};

  auto _ = all_t{};
  for (auto const &w : wmesh)
    chi_wr[w, _] = triqs::gfs::fourier(chi_wk[w, _]);

  return chi_wr;
}
  */
  
} // namespace triqs_tprf
