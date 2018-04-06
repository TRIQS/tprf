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

#include "linalg.hpp"
#include "lattice.hpp"

#include <triqs/arrays/linalg/eigenelements.hpp>

#include <triqs/clef.hpp>
using namespace triqs::clef;

namespace {
placeholder<0> iw;
placeholder<1> inu;
placeholder<2> k;
placeholder<3> r;

placeholder<4> a;
placeholder<5> b;
placeholder<6> c;
placeholder<7> d;

placeholder<8> inup;
placeholder<9> tau;

} // namespace

#include <triqs/arrays/linalg/det_and_inverse.hpp>
using triqs::arrays::inverse;

namespace tprf {

// ----------------------------------------------------
// g

gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh) {
  gk_iw_t g0k = make_gf<gk_iw_t::mesh_t::var_t>({mesh, ek.mesh()}, ek.target());

  g0k(inu, k)(a, b) << kronecker(a, b) * (inu + mu) - ek(k)(a, b);

  for (auto const &kidx : std::get<1>(g0k.mesh())) {
    auto _ = var_t{};
    g0k[_, kidx] = inverse(g0k[_, kidx]);
  }

  return g0k;
}

gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma) {

  gk_iw_t gk =
      make_gf<gk_iw_t::mesh_t::var_t>({sigma.mesh(), ek.mesh()}, ek.target());

  gk(inu, k)(a, b) << kronecker(a, b) * (inu + mu) - ek(k)(a, b) -
                          sigma(inu)(a, b);

  for (auto const &kidx : std::get<1>(gk.mesh())) {
    auto _ = var_t{};
    gk[_, kidx] = inverse(gk[_, kidx]);
  }

  // gk = inverse(gk);  // does not work, see TRIQS issue #463
  return gk;
}

gr_iw_t gr_from_gk(gk_iw_vt gk, gf_mesh<cyclic_lattice> lmesh) {
  int nk = std::get<1>(gk.mesh()).get_dimensions()[0];

  gr_iw_t gr = make_gf<gr_iw_t::mesh_t::var_t>(
      {std::get<0>(gk.mesh()), lmesh}, gk.target());

  auto _ = var_t{};
  for (auto const &nidx : std::get<0>(gr.mesh()))
    gr[nidx, _] = inverse_fourier(gk[nidx, _]);

  return gr;
}

gk_iw_t gk_from_gr(gr_iw_vt gr, brillouin_zone bz) {
  int nk = std::get<1>(gr.mesh()).get_dimensions()[0];

  gk_iw_t gk = make_gf<gk_iw_t::mesh_t::var_t>(
      {std::get<0>(gr.mesh()), {bz, nk}}, gr.target());

  auto _ = var_t{};
  for (auto const &nidx : std::get<0>(gk.mesh()))
    gk[nidx,_] = fourier(gr[nidx,_]);

  return gk;
}

// ----------------------------------------------------
// Transformations: Matsubara frequency <-> imaginary time

/** Set the two lowest order tail coefficients
    
    for a single-particle Green's function.
    using the 1/w first order coefficient
    and fitting the second order coefficient
    from the value at the lowest Matsubara frequency.
 */
  
g_iw_vt set_gf_tail(g_iw_vt gw, double s_1 = 1.) {

  auto s = gw.singularity();
  
  s.zero();
  
  s(1) = s_1; // set order 1 to the unit matrix

  // get approx of 2nd order from lowest G(w) value
  
  for (auto const &w : gw.mesh()) {
    s(2) = w*w*(gw[w] - s(1)/w);
    break;
  }

  return gw;
}
  
gr_tau_t grt_from_grw(gr_iw_vt grw) {

  auto wmesh = std::get<0>(grw.mesh());
  auto rmesh = std::get<1>(grw.mesh());
  
  double beta = wmesh.domain().beta;
  
  int nw = wmesh.last_index() + 1;
  int ntau = 4 * nw;

  gr_tau_t grt = make_gf<gr_tau_t::mesh_t::var_t>(
    {{beta, Fermion, ntau}, rmesh}, grw.target());

  auto _ = var_t{};
  for (auto const &r : rmesh) {

    auto gw =  grw[_, r];

    //std::cout << "r = " << r << "\n";
    //std::cout << "ridx = " << r.linear_index() << "\n";

    double s_1 = 0.;
    if(r.linear_index() == 0) s_1 = 1.; // only r=0 has a fermi jump and a 1/i\omega_n tail
    //std::cout << "s_1 = " << s_1 << "\n";
      
    grt[_, r] = inverse_fourier(set_gf_tail(grw[_, r], s_1)); // Using on the fly tail coeff FIXME

    //grt[_, r] = inverse_fourier(set_gf_tail(grw[_, r])); // Using on the fly tail coeff FIXME
  }

  return grt;
}

// ----------------------------------------------------
// chi00 bubble in analytic form

double fermi(double e) { return 1./(exp(e) + 1.); }
  
chi_wk_t chi00_wk_from_ek(gf<brillouin_zone, matrix_valued> ek_in, int nw, double beta, double mu) {

  auto kmesh = ek_in.mesh();
  int nb = ek_in.target().shape()[0];

  chi_wk_t chi{{{beta, Boson, nw}, kmesh}, {nb, nb, nb, nb}};
  for( auto const & [w, k] : chi.mesh() ) chi[w, k] = 0.;

  auto wmesh = std::get<0>(chi.mesh());
  
  for( auto const &k : kmesh) {

    std::cout << "k = " << k << "\n";

    matrix<std::complex<double>> ek_mat(ek_in[k] - mu);
    
    auto eig_k = linalg::eigenelements(ek_mat);
    auto ek = eig_k.first;
    auto Uk = eig_k.second;

    for( auto const &w : wmesh ) {
      for( auto const &q : kmesh) {
      
       matrix<std::complex<double>> ekq_mat(ek_in(k + q) - mu);
       auto eig_kq = linalg::eigenelements(ekq_mat);
       auto ekq = eig_kq.first;
       auto Ukq = eig_kq.second;

       for( int i : range(nb) )
	 for( int j : range(nb) ) {

	   double de = ekq(j) - ek(i);
	   double dn = fermi(ek(i) * beta) - fermi(ekq(j) * beta);
	   double tol = 1e-10;
	   
	   if( abs(std::complex<double>(w)) < tol && abs(de) < tol ) { // w=0, de=0, 2nd order pole

	     // -- analytic first derivative of the fermi distribution function
	     // -- evaluated at ek(i)
	     
	     double cosh_be = cosh(0.5 * beta * ek(i));
	     double dn_w0 = beta / (4. * cosh_be * cosh_be);

	     chi[w, q](a, b, c, d) << chi[w, q](a, b, c, d) +
	       dagger(Uk)(a, i) * Uk(i, b) * dagger(Ukq)(c, j) * Ukq(j, d) * dn_w0;

	   } else { 

	     chi[w, q](a, b, c, d) << chi[w, q](a, b, c, d) +
	       dagger(Uk)(a, i) * Uk(i, b) * dagger(Ukq)(c, j) * Ukq(j, d) * dn / ( w + de );

	   }
	 }
      }
    }
  }

  chi /= nb * kmesh.size();

  /*
  for( auto const &q : kmesh) {
    std::cout << "q = " << q << "\n";
    std::cout << chi[Idx(0), q] << "\n";
  }
  */
  
  
  return chi;
}

// ----------------------------------------------------
// chi0 bubble in imaginary time

chi_tr_t chi0_tr_from_grt_PH(gr_tau_vt grt) {

  auto tmesh = std::get<0>(grt.mesh());
  auto rmesh = std::get<1>(grt.mesh());

  int nb = grt.target().shape()[0];
  int ntau = tmesh.size();
  double beta = tmesh.domain().beta;

  chi_tr_t chi0_tr{{{beta, Boson, ntau}, rmesh}, {nb, nb, nb, nb}};

  // -- This does not work on the boundaries!! The eval wraps to the other regime!
  // -- gt(beta) == gt(beta + 0^+)    
  //chi0_tr(tau, r)(a, b, c, d) << grt(tau, r)(d, a) * grt(-tau, -r)(b, c);

  for (auto const &r : rmesh) {
    for (auto const &t : tmesh) {

      // -- This does not work on the boundaries!! The eval wraps to the other regime!
      // -- gt(beta) == gt(beta + 0^+)    
      //chi0_tr[t, r](a, b, c, d) << grt(t, r)(d, a) * grt(-t, -r)(b, c);

      // -- Ugly hack to evaluate within the range of t in [0, beta] and -t in [-beta, 0] respectively
      
      double t_p = float(+t);
      double t_m = float(-t);

      double eps = 1e-9;

      // -- Shift any point on the boundary by eps inside the boundary... :P
      
      if(abs(t_p) < eps) t_p = eps;
      if(abs(t_p - beta) < eps) t_p = + beta - eps;

      if(abs(t_m) < eps) t_m = - eps;
      if(abs(t_m + beta) < eps) t_m = - beta + eps;
      
      chi0_tr[t, r](a, b, c, d) << - grt(t_p, r)(d, a) * grt(t_m, -r)(b, c);

    }
  }

  return chi0_tr;
}

chi_wr_t chi_wr_from_chi_tr(chi_tr_vt chi_tr) {

  int nb = chi_tr.target().shape()[0];
  //auto target = chi_tr.target();
  
  auto tmesh = std::get<0>(chi_tr.mesh());
  auto rmesh = std::get<1>(chi_tr.mesh());

  int ntau = tmesh.size();
  int nw = ntau / 4;
  double beta = tmesh.domain().beta;
  
  //chi_wr_t chi_wr{{{beta, Boson, nw}, rmesh}, target};
  chi_wr_t chi_wr{{{beta, Boson, nw}, rmesh}, {nb, nb, nb, nb}};

  auto _ = var_t{};
  for (auto const &r : rmesh) {
    std::cout << "r = " << r << "\n";
    auto chi_t = chi_tr[_, r];

    auto s = chi_t.singularity();
    s.zero();
    //s(1) = 1.;

    chi_wr[_, r] = fourier(chi_t);
    //chi_wr[_, r] = fourier(chi_tr[_, r]);
  }

  return chi_wr;
}

chi_wk_t chi_wk_from_chi_wr(chi_wr_vt chi_wr, gf_mesh<brillouin_zone> kmesh) {

  //auto target = chi_wr.target();
  int nb = chi_wr.target().shape()[0];
  
  auto wmesh = std::get<0>(chi_wr.mesh());
  auto rmesh = std::get<1>(chi_wr.mesh());

  chi_wk_t chi_wk{{wmesh, kmesh}, {nb, nb, nb, nb}};

  auto _ = var_t{};
  for (auto const &w : wmesh)
    chi_wk[w, _] = fourier(chi_wr[w, _]);

  return chi_wk;
}
  
// ----------------------------------------------------
// chi0 bubble in Matsubara frequency

chi0r_t chi0r_from_gr_PH(int nw, int nnu, gr_iw_vt gr) {

  int nb = gr.target().shape()[0];
  auto clmesh = std::get<1>(gr.mesh());
  double beta = std::get<0>(gr.mesh()).domain().beta;

  chi0r_t chi0r{{{beta, Boson, nw}, {beta, Fermion, nnu}, clmesh},
                {nb, nb, nb, nb}};

  chi0r(iw, inu, r)(a, b, c, d) << - beta * gr(inu, r)(d, a) * gr(inu + iw, -r)(b, c);

  return chi0r;
}

chi0r_t chi0r_from_chi0q(chi0q_vt chi0q, gf_mesh<cyclic_lattice> clmesh) {

  auto mb = std::get<0>(chi0q.mesh());
  auto mf = std::get<1>(chi0q.mesh());
  auto bzmesh = std::get<2>(chi0q.mesh());

  auto bl = bzmesh.domain().lattice();
  auto clmesh_ref = gf_mesh<cyclic_lattice>(bl, bzmesh.periodization_matrix());
  
  auto chi0r =
      make_gf<chi0r_t::mesh_t::var_t>({mb, mf, clmesh}, chi0q.target());

  for (auto const &widx : std::get<0>(chi0r.mesh())) {
    for (auto const &nidx : std::get<1>(chi0r.mesh())) {
      auto _ = var_t{};
      chi0r[widx, nidx, _] = inverse_fourier(chi0q[widx, nidx, _]);
    }
  }
  return chi0r;
}

chi0q_t chi0q_from_chi0r(chi0r_vt chi0r, gf_mesh<brillouin_zone> bzmesh) {
  auto mb = std::get<0>(chi0r.mesh());
  auto mf = std::get<1>(chi0r.mesh());

  auto chi0q =
      make_gf<chi0q_t::mesh_t::var_t>({mb, mf, bzmesh}, chi0r.target());

  for (auto const &widx : std::get<0>(chi0q.mesh())) {
    for (auto const &nidx : std::get<1>(chi0q.mesh())) {
      auto _ = var_t{};
      chi0q[widx, nidx, _] = fourier(chi0r[widx, nidx, _]);
    }
  }
  return chi0q;
}

gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>>
chi0q_sum_nu(chi0q_t chi0q) {

  auto mesh = std::get<1>(chi0q.mesh());
  auto chi0q_w = make_gf<cartesian_product<imfreq, brillouin_zone>>(
      {std::get<0>(chi0q.mesh()), std::get<2>(chi0q.mesh())}, chi0q.target());

  double beta = mesh.domain().beta;
  chi0q_w(iw, k) << sum(chi0q(iw, inu, k), inu = mesh) / (beta*beta);
  return chi0q_w;
}

gf<imfreq, tensor_valued<4>> chi0q_sum_nu_q(chi0q_t chi0q) {

  auto mesh_b = std::get<0>(chi0q.mesh());
  auto mesh_f = std::get<1>(chi0q.mesh());
  auto mesh_k = std::get<2>(chi0q.mesh());
  
  auto chi0_w = make_gf<imfreq>(mesh_b, chi0q.target());

  for(auto const &[w, n, k] : chi0q.mesh())
    chi0_w[w] += chi0q[w, n, k];

  double nk = mesh_k.size();
  double beta = mesh_f.domain().beta;
  chi0_w = chi0_w / nk / (beta*beta);
  
  return chi0_w;
}

  
// ----------------------------------------------------
// chi

chiq_t chiq_from_chi0q_and_gamma_PH(chi0q_vt chi0q, g2_iw_vt gamma_ph) {

  auto _ = var_t{};
  
  auto mb = std::get<0>(chi0q.mesh());
  auto mf = std::get<1>(chi0q.mesh());
  auto mbz = std::get<2>(chi0q.mesh());

  auto chiq =
    make_gf<chiq_t::mesh_t::var_t>({mbz, mb, mf, mf}, chi0q.target());

  for (auto const &k : mbz) {

    // -- Construct matrix version of chi0q_k

    // -- If we could make this a 1,1,1 g2_iw_t function and do the PH inverse
    // -- only in the target space we would save one global inverse! 
    
    auto chi0q_k =
      make_gf<g2_iw_t::mesh_t::var_t>({mb, mf, mf}, chi0q.target());

    for (auto const &w : mb) {
      for (auto const &n : mf) {
	chi0q_k[w, n, n] = chi0q[w, n, k];
      }
    }

    g2_iw_t chiq_inv_k = inverse<Channel_t::PH>(chi0q_k) - gamma_ph;

    chiq[k, _, _, _] = inverse<Channel_t::PH>(chiq_inv_k);
  }
  
  return chiq;
}

gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu(chiq_t chiq) {

  auto mesh_k = std::get<0>(chiq.mesh());
  auto mesh_b = std::get<1>(chiq.mesh());
  auto mesh_f = std::get<2>(chiq.mesh());
  auto chiq_w = make_gf<cartesian_product<brillouin_zone, imfreq>>({mesh_k, mesh_b}, chiq.target());

  // Does not compile due to treatment of the tail (singularity)
  //chiq_w(k, iw) << sum(chiq(k, iw, inu, inup), inu=mesh, inup=mesh);

  for(auto const &[k, w, n1, n2] : chiq.mesh())
    chiq_w[k, w] += chiq[k, w, n1, n2];

  double beta = mesh_f.domain().beta;
  chiq_w = chiq_w / (beta*beta);
  
  return chiq_w;
}

gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(chiq_t chiq) {

  auto mesh_k = std::get<0>(chiq.mesh());
  auto mesh_b = std::get<1>(chiq.mesh());
  auto mesh_f = std::get<2>(chiq.mesh());
  auto chi_w = make_gf<imfreq>(mesh_b, chiq.target());

  for(auto const &[k, w, n1, n2] : chiq.mesh())
    chi_w[w] += chiq[k, w, n1, n2];

  double nk = mesh_k.size();
  double beta = mesh_f.domain().beta;
  chi_w = chi_w / nk / (beta*beta);
  
  return chi_w;
}

} // namespace tprf
