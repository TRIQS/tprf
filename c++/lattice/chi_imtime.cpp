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
#include "chi_imtime.hpp"

namespace tprf {

// ----------------------------------------------------
// chi0 bubble in imaginary time

chi_tr_t chi0_tr_from_grt_PH(gr_tau_vt grt) {

  auto _ = var_t{};

  auto tmesh = std::get<0>(grt.mesh());
  auto rmesh = std::get<1>(grt.mesh());

  int nb = grt.target().shape()[0];
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_tr_t chi0_tr{{{beta, Boson, ntau}, rmesh}, {nb, nb, nb, nb}};

  auto g_target = grt.target();
  auto chi_target = chi0_tr.target();
  
  // -- This does not work on the boundaries!! The eval wraps to the other
  // regime!
  // -- gt(beta) == gt(beta + 0^+)
  // chi0_tr(tau, r)(a, b, c, d) << grt(tau, r)(d, a) * grt(-tau, -r)(b, c);

  //for (auto const &r : rmesh) {

#pragma omp parallel for 
  for (int idx = 0; idx < rmesh.size(); idx++) {
    auto iter = rmesh.begin(); iter += idx; auto r = *iter;

    auto chi0_t = make_gf<imtime>({beta, Boson, ntau}, chi_target);
    auto g_pr_t = make_gf<imtime>(tmesh, g_target);
    auto g_mr_t = make_gf<imtime>(tmesh, g_target);

#pragma omp critical
    {
      g_pr_t = grt[_, r];
      g_mr_t = grt[_, -r];
    }
    
    for (auto const &t : tmesh) {

      // -- This does not work on the boundaries!! The eval wraps to the other
      // regime!
      // -- gt(beta) == gt(beta + 0^+)
      //chi0_tr[t, r](a, b, c, d) << grt(t, r)(d, a) * grt(-t, -r)(b, c);

      chi0_t[t](a, b, c, d) << g_pr_t(t)(d, a) * g_mr_t(beta - t)(b, c);

      /*
      // -- Ugly hack to evaluate within the range of t in [0, beta] and -t in
      // [-beta, 0] respectively

      double t_p = double(+t);
      double t_m = double(-t);

      double eps = 1e-9;

      // -- Shift any point on the boundary by eps inside the boundary... :P

      if (abs(t_p) < eps)
        t_p = eps;
      if (abs(t_p - beta) < eps)
        t_p = +beta - eps;

      if (abs(t_m) < eps)
        t_m = -eps;
      if (abs(t_m + beta) < eps)
        t_m = -beta + eps;

      //chi0_tr[t, r](a, b, c, d) << -grt(t_p, r)(d, a) * grt(t_m, -r)(b, c);
      
      chi0_t[t](a, b, c, d) << -g_pr_t(t_p)(d, a) * g_mr_t(t_m)(b, c);
      */
    }

#pragma omp critical
    chi0_tr[_, r] = chi0_t;
    
  }

  return chi0_tr;
}

// -- optimized version for w=0
chi_wr_t chi0_w0r_from_grt_PH(gr_tau_vt grt) {

  auto _ = var_t{};

  auto tmesh = std::get<0>(grt.mesh());
  auto rmesh = std::get<1>(grt.mesh());

  int nw = 1;
  int nb = grt.target().shape()[0];
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_wr_t chi0_wr{{{beta, Boson, nw}, rmesh}, {nb, nb, nb, nb}};

  auto g_target = grt.target();
  auto chi_target = chi0_wr.target();
  
  //for (auto const &r : rmesh) {

#pragma omp parallel for 
  for (int idx = 0; idx < rmesh.size(); idx++) {
    auto iter = rmesh.begin(); iter += idx; auto r = *iter;

    auto chi0_t = make_gf<imtime>({beta, Boson, ntau}, chi_target);
    auto g_pr_t = make_gf<imtime>(tmesh, g_target);
    auto g_mr_t = make_gf<imtime>(tmesh, g_target);

#pragma omp critical
    {
      g_pr_t = grt[_, r];
      g_mr_t = grt[_, -r];
    }
    
    for (auto const &t : tmesh) {
      chi0_t[t](a, b, c, d) << g_pr_t(t)(d, a) * g_mr_t(beta - t)(b, c);
    }

    auto int_chi0 = chi_trapz_tau(chi0_t);
    
#pragma omp critical
    chi0_wr[0, r] = int_chi0;

  }

  return chi0_wr;
}  

chi_t_t::zero_t chi_trapz_tau(chi_t_vt chi_t) {

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
chi_wr_t chi_w0r_from_chi_tr(chi_tr_vt chi_tr) {

  int nb = chi_tr.target().shape()[0];

  auto tmesh = std::get<0>(chi_tr.mesh());
  auto rmesh = std::get<1>(chi_tr.mesh());

  int nw = 1;
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_wr_t chi_wr{{{beta, Boson, nw}, rmesh}, {nb, nb, nb, nb}};

  //for (auto const &r : rmesh) {

#pragma omp parallel for 
  for (int idx = 0; idx < rmesh.size(); idx++) {
    auto iter = rmesh.begin(); iter += idx; auto r = *iter;

    /*
    chi_wr[0, r] *= 0.;
    
    // -- Trapetzoidal integration by hand
    for (auto const &t : tmesh)
      chi_wr[0, r] += chi_tr[t, r];

    auto boundary_iter = tmesh.begin();
    chi_wr[0, r] -= 0.5 * chi_tr[*boundary_iter, r];
    boundary_iter += ntau - 1;
    chi_wr[0, r] -= 0.5 * chi_tr[*boundary_iter, r];

    chi_wr[0, r] *= beta / (ntau-1);
    */

    auto _ = var_t{};
    chi_wr[0, r] = chi_trapz_tau(chi_tr[_, r]);
  }

  return chi_wr;
}

chi_wr_t chi_wr_from_chi_tr(chi_tr_vt chi_tr, int nw) {

  int nb = chi_tr.target().shape()[0];
  // auto target = chi_tr.target();

  auto tmesh = std::get<0>(chi_tr.mesh());
  auto rmesh = std::get<1>(chi_tr.mesh());

  int ntau = tmesh.size();
  //int nw = ntau / 4;
  double beta = tmesh.domain().beta;

  // chi_wr_t chi_wr{{{beta, Boson, nw}, rmesh}, target};
  chi_wr_t chi_wr{{{beta, Boson, nw}, rmesh}, {nb, nb, nb, nb}};

  for (auto const &r : rmesh) {
  /*
#pragma omp parallel for 
  for (int idx = 0; idx < rmesh.size(); idx++) {
    auto iter = rmesh.begin(); iter += idx; auto r = *iter;
  */

    //std::cout << "r = " << r << "\n";

    /*
    auto _ = var_t{};
    auto chi_t = chi_tr[_, r];

    auto s = chi_t.singularity();
    s.zero();
    // s(1) = 1.;

    chi_wr[_, r] = fourier(chi_t);
    */

    auto _ = var_t{};
    chi_wr[_, r] = fourier(chi_tr[_, r]);
  }

  return chi_wr;
}
  
chi_wk_t chi_wk_from_chi_wr(chi_wr_vt chi_wr) {

  // auto target = chi_wr.target();
  int nb = chi_wr.target().shape()[0];

  auto wmesh = std::get<0>(chi_wr.mesh());
  auto rmesh = std::get<1>(chi_wr.mesh());

  auto kmesh = gf_mesh<brillouin_zone>{brillouin_zone{rmesh.domain()}, rmesh.periodization_matrix};
  
  chi_wk_t chi_wk{{wmesh, kmesh}, {nb, nb, nb, nb}};

  auto _ = var_t{};
  for (auto const &w : wmesh)
    chi_wk[w, _] = fourier(chi_wr[w, _]);

  return chi_wk;
}

} // namespace tprf
