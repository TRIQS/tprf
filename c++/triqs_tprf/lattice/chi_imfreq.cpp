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

#include <triqs/utility/timer.hpp>
#include <triqs/utility/timestamp.hpp>
#include <omp.h>

#include "../fourier/fourier.hpp"
#include "../linalg.hpp"
#include "../mpi.hpp"

#include "chi_imfreq.hpp"
#include "common.hpp"

namespace triqs_tprf {

  namespace {
    using fourier::_fourier_plan;
    using fourier::_fourier_with_plan;
  } // namespace

  // ----------------------------------------------------
  // chi0 bubble in Matsubara frequency

  /*
chi0r_t chi0r_from_gr_PH(int nw, int nnu, gr_iw_vt gr) {

int nb = gr.target().shape()[0];
auto clmesh = std::get<1>(gr.mesh());
double beta = std::get<0>(gr.mesh()).domain().beta;

chi0r_t chi0r{{{beta, Boson, nw}, {beta, Fermion, nnu}, clmesh},
              {nb, nb, nb, nb}};

chi0r(iw, inu, r)(a, b, c, d)
    << -beta * gr(inu, r)(d, a) * gr(inu + iw, -r)(b, c);

return chi0r;
}
*/

  chi_wnr_t chi0r_from_gr_PH(int nw, int nn, g_wr_cvt g_nr) {

    mpi::communicator comm;

    triqs::utility::timer t_alloc, t_calc, t_mpi_all_reduce;

    auto _ = all_t{};

    int nb     = g_nr.target().shape()[0];
    auto rmesh = std::get<1>(g_nr.mesh());

    double beta = std::get<0>(g_nr.mesh()).domain().beta;

    auto wmesh = gf_mesh<imfreq>{beta, Boson, nw};
    auto nmesh = gf_mesh<imfreq>{beta, Fermion, nn};

    t_alloc.start();
    chi0r_t chi0_wnr{{wmesh, nmesh, rmesh}, {nb, nb, nb, nb}};
    chi0_wnr *= 0.;
    t_alloc.stop();

    auto chi_target = chi0_wnr.target();
    auto g_target   = g_nr.target();

    auto arr = mpi_view(rmesh);

    std::cout << "rank " << comm.rank() << " has arr.size() = " << arr.size() << " of " << rmesh.size() << std::endl;

    t_calc.start();
#pragma omp parallel for
    for (int idx = 0; idx < arr.size(); idx++) {
      auto &r = arr(idx);

      auto chi0_wn = make_gf<cartesian_product<imfreq, imfreq>>({wmesh, nmesh}, chi_target);
      auto g_pr_n  = make_gf<imfreq>(nmesh, g_target);
      auto g_mr_n  = make_gf<imfreq>(nmesh, g_target);

#pragma omp critical
      {
        g_pr_n = g_nr[_, r];
        g_mr_n = g_nr[_, -r];
      }

      for (auto const &w : wmesh)
        for (auto const &n : nmesh) chi0_wn[w, n](a, b, c, d) << -beta * g_pr_n(n)(d, a) * g_mr_n(n + w)(b, c);

#pragma omp critical
      chi0_wnr[_, _, r] = chi0_wn;

      // chi0r(iw, inu, r)(a, b, c, d) << -beta * gr(inu, r)(d, a) * gr(inu + iw,
      // -r)(b, c);
    }
    t_calc.stop();

    t_mpi_all_reduce.start();

    // This (surprisingly gives incorrect results, for sufficiently large arguments...
    //chi0_wnr = mpi::all_reduce(chi0_wnr);

    /* // Does not give correct result either...
  for( auto const & w : wmesh )
    chi0_wnr[w, _, _] = mpi::all_reduce(chi0_wnr[w, _, _]);
  */

    for (auto const &w : wmesh)
      for (auto const &n : nmesh) chi0_wnr[w, n, _] = mpi::all_reduce(chi0_wnr[w, n, _]);

    t_mpi_all_reduce.stop();

    std::cout << "--> chi0r_from_gr_PH: alloc " << double(t_alloc) << " s, calc " << double(t_calc) << " s, mpi::all_reduce "
              << double(t_mpi_all_reduce) << " s." << std::endl;

    return chi0_wnr;
  }

  // ---

  chi_wnr_t chi0r_from_gr_PH_nompi(int nw, int nn, g_wr_cvt g_nr) {

    triqs::utility::timer t_alloc, t_calc, t_mpi_all_reduce;

    auto _ = all_t{};

    int nb     = g_nr.target().shape()[0];
    auto rmesh = std::get<1>(g_nr.mesh());

    double beta = std::get<0>(g_nr.mesh()).domain().beta;

    auto wmesh = gf_mesh<imfreq>{beta, Boson, nw};
    auto nmesh = gf_mesh<imfreq>{beta, Fermion, nn};

    t_alloc.start();
    chi0r_t chi0_wnr{{wmesh, nmesh, rmesh}, {nb, nb, nb, nb}};
    chi0_wnr *= 0.;
    t_alloc.stop();

    auto chi_target = chi0_wnr.target();
    auto g_target   = g_nr.target();

    t_calc.start();

#pragma omp parallel for
    for (int idx = 0; idx < rmesh.size(); idx++) {
      auto iter = rmesh.begin();
      iter += idx;
      auto r = *iter;

      auto chi0_wn = make_gf<cartesian_product<imfreq, imfreq>>({wmesh, nmesh}, chi_target);
      auto g_pr_n  = make_gf<imfreq>(nmesh, g_target);
      auto g_mr_n  = make_gf<imfreq>(nmesh, g_target);

#pragma omp critical
      {
        g_pr_n = g_nr[_, r];
        g_mr_n = g_nr[_, -r];
      }

      for (auto const &w : wmesh)
        for (auto const &n : nmesh) chi0_wn[w, n](a, b, c, d) << -beta * g_pr_n(n)(d, a) * g_mr_n(n + w)(b, c);

#pragma omp critical
      chi0_wnr[_, _, r] = chi0_wn;

      // chi0r(iw, inu, r)(a, b, c, d) << -beta * gr(inu, r)(d, a) * gr(inu + iw,
      // -r)(b, c);
    }
    t_calc.stop();

    t_mpi_all_reduce.start();
    //chi0_wnr = mpi::all_reduce(chi0_wnr);
    t_mpi_all_reduce.stop();

    std::cout << "--> chi0r_from_gr_PH: alloc " << double(t_alloc) << " s, calc " << double(t_calc) << " s, mpi::all_reduce "
              << double(t_mpi_all_reduce) << " s." << std::endl;

    return chi0_wnr;
  }

  // ----------------------------------------------------

  // Helper function compiting chi0 for fixed bosonic frequency w and momentum q.

  CPP2PY_IGNORE
  gf<imfreq, tensor_valued<4>> chi0_n_from_g_wk_PH(mesh_point<gf_mesh<imfreq>> w, mesh_point<cluster_mesh> q, gf_mesh<imfreq> fmesh, g_wk_cvt g_wk) {

    int nb                    = g_wk.target().shape()[0];
    auto [fmesh_large, kmesh] = g_wk.mesh();

    double beta = fmesh.domain().beta;

    gf<imfreq, tensor_valued<4>> chi0_n{fmesh, {nb, nb, nb, nb}};

    // 100x times slower
    // chi0_n(inu)(a, b, c, d) << -beta/kmesh.size() * sum(g_wk(inu, k)(d, a) *
    // g_wk(inu + w, k - q)(b, c), k=kmesh);

    for (auto const &n : fmesh) {
      for (auto const &k : kmesh) {
        auto g_da = g_wk[n, k];
        auto g_bc = g_wk[n + w, k - q];
        for (auto a : range(nb))
          for (auto b : range(nb))
            for (auto c : range(nb))
              for (auto d : range(nb)) chi0_n[n](a, b, c, d) -= g_da(d, a) * g_bc(b, c);
      }
    }

    chi0_n *= beta / kmesh.size();

    return chi0_n;
  }

  // ----------------------------------------------------

  // Helper function compiting chi0 for fixed bosonic frequency w and momentum q.
  // using the self energy and the dispersion (instead of the greens function)

  CPP2PY_IGNORE
  gf<imfreq, tensor_valued<4>> chi0_n_from_e_k_sigma_w_PH(mesh_point<gf_mesh<imfreq>> w, mesh_point<cluster_mesh> q, gf_mesh<imfreq> fmesh, double mu,
                                                          e_k_cvt e_k, g_w_cvt sigma_w) {

    int nb     = e_k.target().shape()[0];
    auto kmesh = e_k.mesh();

    auto fmesh_large = sigma_w.mesh();

    double beta = fmesh.domain().beta;
    auto I      = make_unit_matrix<ek_vt::scalar_t>(e_k.target_shape()[0]);

    gf<imfreq, tensor_valued<4>> chi0_n{fmesh, {nb, nb, nb, nb}};

    for (auto const &k : kmesh) {
      for (auto const &n : fmesh) {

        auto g_da = inverse((n + mu) * I - e_k[k] - sigma_w[matsubara_freq(n)]);
        auto g_bc = inverse((n + mu) * I - e_k[k - q] - sigma_w[n + w]);

        for (auto a : range(nb))
          for (auto b : range(nb))
            for (auto c : range(nb))
              for (auto d : range(nb)) chi0_n[n](a, b, c, d) -= g_da(d, a) * g_bc(b, c);
      }
    }

    chi0_n *= beta / kmesh.size();

    return chi0_n;
  }

  // ----------------------------------------------------

  chi_wnk_t chi0q_from_g_wk_PH(int nw, int nn, g_wk_cvt g_wk) {

    auto [fmesh_large, kmesh] = g_wk.mesh();

    int nb      = g_wk.target().shape()[0];
    double beta = std::get<0>(g_wk.mesh()).domain().beta;

    gf_mesh<imfreq> bmesh{beta, Boson, nw};
    gf_mesh<imfreq> fmesh{beta, Fermion, nn};

    assert(fmesh.size() < fmesh_large.size());

    chi0q_t chi0_wnk({bmesh, fmesh, kmesh}, {nb, nb, nb, nb});

    auto _ = all_t{};
    for (auto const &[w, q] : gf_mesh{bmesh, kmesh}) { chi0_wnk[w, _, q] = chi0_n_from_g_wk_PH(w, q, fmesh, g_wk); }

    return chi0_wnk;
  }

  // ----------------------------------------------------
  // momentum <-> realspace transforms

  /*
chi0r_t chi0r_from_chi0q(chi0q_vt chi0_wnk) {

auto [bmesh, fmesh, kmesh] = chi0_wnk.mesh();
auto rmesh = make_adjoint_mesh(kmesh);

auto chi0_wnr =
    make_gf<chi0r_t::mesh_t::var_t>({bmesh, fmesh, rmesh}, chi0_wnk.target());

auto _ = all_t{};
for (auto const &[w, n] : mpi_view(gf_mesh{bmesh, fmesh}))
  chi0_wnr[w, n, _] = fourier(chi0_wnk[w, n, _]);
chi0_wnr = mpi::all_reduce(chi0_wnr);

return chi0_wnr;
}
*/

  chi_wnr_t chi0r_from_chi0q(chi_wnk_cvt chi_wnk) {

    auto _      = all_t{};
    auto target = chi_wnk.target();

    //auto [bmesh, fmesh, kmesh] = chi0_wnk.mesh(); // clang+OpenMP can not handle this...
    auto bmesh = std::get<0>(chi_wnk.mesh());
    auto fmesh = std::get<1>(chi_wnk.mesh());
    auto kmesh = std::get<2>(chi_wnk.mesh());
    auto rmesh = make_adjoint_mesh(kmesh);

    chi_wnr_t chi_wnr({bmesh, fmesh, rmesh}, chi_wnk.target_shape());

    auto w0 = *bmesh.begin();
    auto n0 = *fmesh.begin();
    auto p  = _fourier_plan<0>(gf_const_view(chi_wnk[w0, n0, _]), gf_view(chi_wnr[w0, n0, _]));

    auto arr = mpi_view(gf_mesh{bmesh, fmesh});

#pragma omp parallel for shared(kmesh, rmesh)
    for (int idx = 0; idx < arr.size(); idx++) {
      //auto &[w, n] = arr(idx);
      auto w = std::get<0>(arr(idx));
      auto n = std::get<1>(arr(idx));

      auto chi_r = make_gf<cyclic_lattice>(rmesh, target);
      auto chi_k = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
      chi_k = chi_wnk[w, n, _];

      _fourier_with_plan<0>(gf_const_view(chi_k), gf_view(chi_r), p);

#pragma omp critical
      chi_wnr[w, n, _] = chi_r;
    }

    //chi_wnr = mpi::all_reduce(chi_wnr); // Incorrect results for large args!!

    // Workaround.. :P
    for (auto const &w : bmesh)
      for (auto const &n : fmesh) chi_wnr[w, n, _] = mpi::all_reduce(chi_wnr[w, n, _]);

    return chi_wnr;
  }

  // --
  /*
chi0q_t chi0q_from_chi0r(chi0r_vt chi0_wnr) {

auto [bmesh, fmesh, rmesh] = chi0_wnr.mesh();
auto kmesh = make_adjoint_mesh(rmesh);

auto chi0_wnk =
    make_gf<chi0q_t::mesh_t::var_t>({bmesh, fmesh, kmesh}, chi0_wnr.target());

auto _ = all_t{};
for (auto const &[w, n] : mpi_view(gf_mesh{bmesh, fmesh}))
  chi0_wnk[w, n, _] = triqs::gfs::fourier(chi0_wnr[w, n, _]);
chi0_wnk = mpi::all_reduce(chi0_wnk);

return chi0_wnk;
}
*/

  chi_wnk_t chi0q_from_chi0r(chi_wnr_cvt chi_wnr) {

    triqs::utility::timer t_alloc, t_calc, t_mpi_all_reduce;

    auto _      = all_t{};
    auto target = chi_wnr.target();

    //auto [bmesh, fmesh, rmesh] = chi_wnr.mesh();
    auto bmesh = std::get<0>(chi_wnr.mesh());
    auto fmesh = std::get<1>(chi_wnr.mesh());
    auto rmesh = std::get<2>(chi_wnr.mesh());
    auto kmesh = make_adjoint_mesh(rmesh);

    t_alloc.start();
    chi_wnk_t chi_wnk({bmesh, fmesh, kmesh}, chi_wnr.target_shape());

    auto w0 = *bmesh.begin();
    auto n0 = *fmesh.begin();
    auto p  = _fourier_plan<0>(gf_const_view(chi_wnr[w0, n0, _]), gf_view(chi_wnk[w0, n0, _]));

    auto arr = mpi_view(gf_mesh{bmesh, fmesh});

    t_alloc.stop();
    t_calc.start();

#pragma omp parallel for
    for (int idx = 0; idx < arr.size(); idx++) {
      //auto &[w, n] = arr(idx);
      auto w = std::get<0>(arr(idx));
      auto n = std::get<1>(arr(idx));

      auto chi_r = make_gf<cyclic_lattice>(rmesh, target);
      auto chi_k = make_gf<brillouin_zone>(kmesh, target);

#pragma omp critical
      chi_r = chi_wnr[w, n, _];

      _fourier_with_plan<0>(gf_const_view(chi_r), gf_view(chi_k), p);

#pragma omp critical
      chi_wnk[w, n, _] = chi_k;
    }

    t_calc.stop();
    t_mpi_all_reduce.start();

    //chi_wnk = mpi::all_reduce(chi_wnk); // Incorrect results for large args!!

    // Workaround.. :P
    for (auto const &w : bmesh)
      for (auto const &n : fmesh) chi_wnk[w, n, _] = mpi::all_reduce(chi_wnk[w, n, _]);

    t_mpi_all_reduce.stop();

    std::cout << "--> chi0q_from_chi0r: alloc " << double(t_alloc) << " s, calc " << double(t_calc) << " s, mpi::all_reduce "
              << double(t_mpi_all_reduce) << " s." << std::endl;

    return chi_wnk;
  }

  // ----------------------------------------------------

  chi_wk_t chi0q_sum_nu(chi_wnk_cvt chi_wnk) {

    auto wmesh = std::get<0>(chi_wnk.mesh());
    auto nmesh = std::get<1>(chi_wnk.mesh());
    auto kmesh = std::get<2>(chi_wnk.mesh());

    chi_wk_t chi_wk({wmesh, kmesh}, chi_wnk.target_shape());
    chi_wk *= 0.;

    double beta = wmesh.domain().beta;
    auto arr    = mpi_view(gf_mesh{wmesh, kmesh});

#pragma omp parallel for
    for (int idx = 0; idx < arr.size(); idx++) {
      auto &[w, k] = arr(idx);
      for (auto &n : nmesh) chi_wk[w, k] += chi_wnk[w, n, k];
      chi_wk[w, k] /= beta * beta;
    }

    //chi_wk(iw, k) << sum(chi_wnk(iw, inu, k), inu = mesh) / (beta * beta);

    chi_wk = mpi::all_reduce(chi_wk);
    return chi_wk;
  }

  // ----------------------------------------------------

  chi_wk_t chi0q_sum_nu_tail_corr_PH(chi_wnk_cvt chi_wnk) {

    auto wmesh = std::get<0>(chi_wnk.mesh());
    auto nmesh = std::get<1>(chi_wnk.mesh());
    auto qmesh = std::get<2>(chi_wnk.mesh());

    chi_wk_t chi_wk({wmesh, qmesh}, chi_wnk.target_shape());

    int nb      = chi_wnk.target_shape()[0];
    double beta = wmesh.domain().beta;

    auto chi_target = chi_wnk.target();

    // for (auto const &w : wmesh) {
    //  for (auto const &q : qmesh) {

    auto wq_mesh = gf_mesh{wmesh, qmesh};
    auto arr     = mpi_view(wq_mesh); // FIXME Use library implementation

#pragma omp parallel for
    for (int idx = 0; idx < arr.size(); idx++) {
      //auto &[w, q] = arr(idx);
      auto w = std::get<0>(arr(idx));
      auto q = std::get<1>(arr(idx));

      auto _ = all_t{};

      auto chi = make_gf<imfreq>(nmesh, chi_target);
      array<std::complex<double>, 4> dens(nb, nb, nb, nb);

#pragma omp critical
      chi = chi_wnk[w, _, q];

      for (auto a : range(nb)) {
        for (auto b : range(nb)) {
          for (auto c : range(nb)) {
            for (auto d : range(nb)) {
              auto chi_abcd = slice_target_to_scalar(chi, a, b, c, d);
              // chi0q_w[w, q](a, b, c, d) = density(chi_abcd) / beta;
              dens(a, b, c, d) = density(chi_abcd) / beta;
            }
          }
        }
      }

#pragma omp critical
      chi_wk[w, q] = dens;
    }

    chi_wk = mpi::all_reduce(chi_wk);
    return chi_wk;
  }

  // ----------------------------------------------------

  chi_w_t chi0q_sum_nu_q(chi_wnk_cvt chi_wnk) {

    auto mesh_b = std::get<0>(chi_wnk.mesh());
    auto mesh_f = std::get<1>(chi_wnk.mesh());
    auto mesh_k = std::get<2>(chi_wnk.mesh());

    chi_w_t chi_w(mesh_b, chi_wnk.target_shape());

    for (auto const &[w, n, k] : chi_wnk.mesh()) chi_w[w] += chi_wnk[w, n, k];

    double nk   = mesh_k.size();
    double beta = mesh_f.domain().beta;
    chi_w       = chi_w / nk / (beta * beta);

    return chi_w;
  }

  // ----------------------------------------------------
  // chi

  /*
chiq_t chiq_from_chi0q_and_gamma_PH(chi0q_vt chi0q, g2_iw_vt gamma_ph) {

  auto _ = all_t{};

  auto mb = std::get<0>(chi0q.mesh());
  auto mf = std::get<1>(chi0q.mesh());
  auto mbz = std::get<2>(chi0q.mesh());

  auto chiq = make_gf<chiq_t::mesh_t::var_t>({mbz, mb, mf, mf}, chi0q.target());

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
*/

  // ----------------------------------------------------

  chi_kwnn_t chiq_from_chi0q_and_gamma_PH(chi_wnk_cvt chi0_wnk, chi_wnn_cvt gamma_ph_wnn) {

    auto _ = all_t{};

    auto mb  = std::get<0>(chi0_wnk.mesh());
    auto mf  = std::get<1>(chi0_wnk.mesh());
    auto mbz = std::get<2>(chi0_wnk.mesh());

    chi_kwnn_t chi_kwnn({mbz, mb, mf, mf}, chi0_wnk.target_shape());

    // for (auto const &k : mbz) {

#pragma omp parallel for
    for (int idx = 0; idx < mbz.size(); idx++) {
      auto iter = mbz.begin();
      iter += idx;
      auto k = *iter;

      auto chi0 = make_gf<g2_nn_t::mesh_t::var_t>({mf, mf}, chi0_wnk.target());
      auto I    = identity<Channel_t::PH>(chi0);

      for (auto const &w : mb) {

        chi0 *= 0.;
        for (auto const &n : mf) { chi0[n, n] = chi0_wnk[w, n, k]; }

        // this step could be optimized, using the diagonality of chi0 and I
        chi_nn_t denom = I - product<Channel_t::PH>(chi0, gamma_ph_wnn[w, _, _]);

        // also the last product here
        chi_nn_t chi = product<Channel_t::PH>(inverse<Channel_t::PH>(denom), chi0);

#pragma omp critical
        chi_kwnn[k, w, _, _] = chi;
      }
    }

    return chi_kwnn;
  }

  chi_kwnn_t f_from_chi0q_and_gamma_PH(chi_wnk_cvt chi0_wnk, chi_wnn_cvt gamma_ph_wnn) {

    auto _ = all_t{};

    mpi::communicator c;

    auto target_shape = chi0_wnk.target_shape();
    auto bmesh        = std::get<0>(chi0_wnk.mesh());
    auto fmesh        = std::get<1>(chi0_wnk.mesh());
    auto kmesh        = std::get<2>(chi0_wnk.mesh());

    chi_kwnn_t chi_kwnn({kmesh, bmesh, fmesh, fmesh}, target_shape);

    auto kb_mesh = gf_mesh{kmesh, bmesh};
    auto arr     = mpi_view(kb_mesh); // FIXME Use library implementation

    std::cout << "BSE rank " << c.rank() << " of " << c.size() << " has " << arr.size() << " jobs." << std::endl;

    triqs::utility::timer t;
    t.start();

#pragma omp parallel for
    for (int idx = 0; idx < arr.size(); idx++) {
      //auto &[k, w] = arr(idx);
      auto k = std::get<0>(arr(idx));
      auto w = std::get<1>(arr(idx));

      //triqs::utility::timer t_copy_1, t_chi0_nn, t_bse, t_chi_tr, t_copy_2;

      // ----------------------------------------------------
      //t_copy_1.start();

      chi_w_t chi0_n(fmesh, target_shape);

#pragma omp critical
      chi0_n = chi0_wnk[w, _, k];

      //t_copy_1.stop();
      //std::cout << "BSE: copy_1 " << double(t_copy_1) << " s" << std::endl;
      // ----------------------------------------------------
      //t_chi0_nn.start();

      chi_nn_t chi0_nn({fmesh, fmesh}, target_shape);
      chi0_nn *= 0.;

      for (auto const &n : fmesh) chi0_nn[n, n] = chi0_n[n];

      //t_chi0_nn.stop();
      //std::cout << "BSE: chi0_nn " << double(t_chi0_nn) << " s" << std::endl;
      // ----------------------------------------------------
      //t_bse.start();

      auto I = identity<Channel_t::PH>(chi0_nn);

      chi_nn_t gamma_nn({fmesh, fmesh}, target_shape);

#pragma omp critical
      gamma_nn = gamma_ph_wnn[w, _, _];

      // this step could be optimized, using the diagonality of chi0 and I
      chi_nn_t denom = I - product<Channel_t::PH>(gamma_nn, chi0_nn);

      // also the last product here
      chi_nn_t chi_nn = product<Channel_t::PH>(inverse<Channel_t::PH>(denom), gamma_nn);

      //t_bse.stop();
      //std::cout << "BSE: bse inv " << double(t_bse) << " s" << std::endl;
      // ----------------------------------------------------
      //t_chi_tr.start();

      //t_chi_tr.stop();
      //std::cout << "BSE: chi tr " << double(t_chi_tr) << " s" << std::endl;
      // ----------------------------------------------------
      //t_copy_2.start();

#pragma omp critical
      chi_kwnn[k, w, _, _] = chi_nn;

      //t_copy_2.stop();
      //std::cout << "BSE: copy_2 " << double(t_copy_2) << " s" << std::endl;
      // ----------------------------------------------------

      if (c.rank() == 0 && omp_get_thread_num() == 0) {

        int Nomp = omp_get_num_threads();
        int N    = int(floor(arr.size() / double(Nomp)));

        //double t_left = double(t) * ( N / (idx + 1) - 1. );
        int done_percent = int(floor(100 * double(idx + 1) / N));

        std::cout << "BSE " << triqs::utility::timestamp() << " " << std::setfill(' ') << std::setw(3) << done_percent << "% "
                  << "ETA " << triqs::utility::estimate_time_left(N, idx, t) << " job no " << std::setfill(' ') << std::setw(5) << int(idx)
                  << " of approx " << N << " jobs/thread/rank." << std::endl;
      }
    }

    chi_kwnn = mpi::all_reduce(chi_kwnn);

    t.stop();
    if (c.rank() == 0) std::cout << "BSE TIME: " << double(t) << " s" << std::endl;

    return chi_kwnn;
  }
  // ----------------------------------------------------

  chi_kw_t chiq_sum_nu_from_chi0q_and_gamma_PH(chi_wnk_cvt chi0_wnk, chi_wnn_cvt gamma_ph_wnn) {

    auto _ = all_t{};

    mpi::communicator c;

    auto target_shape = chi0_wnk.target_shape();
    auto bmesh        = std::get<0>(chi0_wnk.mesh());
    auto fmesh        = std::get<1>(chi0_wnk.mesh());
    auto kmesh        = std::get<2>(chi0_wnk.mesh());

    double beta = fmesh.domain().beta;

    chi_kw_t chi_kw({kmesh, bmesh}, target_shape);

    auto kb_mesh = gf_mesh{kmesh, bmesh};
    auto arr     = mpi_view(kb_mesh); // FIXME Use library implementation

    std::cout << "BSE rank " << c.rank() << " of " << c.size() << " has " << arr.size() << " jobs." << std::endl;

    triqs::utility::timer t;
    t.start();

#pragma omp parallel for
    for (int idx = 0; idx < arr.size(); idx++) {
      //auto &[k, w] = arr(idx);
      auto k = std::get<0>(arr(idx));
      auto w = std::get<1>(arr(idx));

      //triqs::utility::timer t_copy_1, t_chi0_nn, t_bse, t_chi_tr, t_copy_2;

      // ----------------------------------------------------
      //t_copy_1.start();

      chi_w_t chi0_n(fmesh, target_shape);

#pragma omp critical
      chi0_n = chi0_wnk[w, _, k];

      //t_copy_1.stop();
      //std::cout << "BSE: copy_1 " << double(t_copy_1) << " s" << std::endl;
      // ----------------------------------------------------
      //t_chi0_nn.start();

      chi_nn_t chi0_nn({fmesh, fmesh}, target_shape);
      chi0_nn *= 0.;

      for (auto const &n : fmesh) chi0_nn[n, n] = chi0_n[n];

      //t_chi0_nn.stop();
      //std::cout << "BSE: chi0_nn " << double(t_chi0_nn) << " s" << std::endl;
      // ----------------------------------------------------
      //t_bse.start();

      auto I = identity<Channel_t::PH>(chi0_nn);

      chi_nn_t gamma_nn({fmesh, fmesh}, target_shape);

#pragma omp critical
      gamma_nn = gamma_ph_wnn[w, _, _];

      // this step could be optimized, using the diagonality of chi0 and I
      chi_nn_t denom = I - product<Channel_t::PH>(chi0_nn, gamma_nn);

      // also the last product here
      chi_nn_t chi = product<Channel_t::PH>(inverse<Channel_t::PH>(denom), chi0_nn);

      //t_bse.stop();
      //std::cout << "BSE: bse inv " << double(t_bse) << " s" << std::endl;
      // ----------------------------------------------------
      //t_chi_tr.start();

      // trace out fermionic frequencies
      array<std::complex<double>, 4> tr_chi(target_shape);

      tr_chi *= 0.0;
      for (auto const &n1 : fmesh)
        for (auto const &n2 : fmesh) tr_chi += chi[n1, n2];

      tr_chi /= beta * beta;

      //t_chi_tr.stop();
      //std::cout << "BSE: chi tr " << double(t_chi_tr) << " s" << std::endl;
      // ----------------------------------------------------
      //t_copy_2.start();

#pragma omp critical
      chi_kw[k, w] = tr_chi;

      //t_copy_2.stop();
      //std::cout << "BSE: copy_2 " << double(t_copy_2) << " s" << std::endl;
      // ----------------------------------------------------

      if (c.rank() == 0 && omp_get_thread_num() == 0) {

        int Nomp = omp_get_num_threads();
        int N    = int(floor(arr.size() / double(Nomp)));

        //double t_left = double(t) * ( N / (idx + 1) - 1. );
        int done_percent = int(floor(100 * double(idx + 1) / N));

        std::cout << "BSE " << triqs::utility::timestamp() << " " << std::setfill(' ') << std::setw(3) << done_percent << "% "
                  << "ETA " << triqs::utility::estimate_time_left(N, idx, t) << " job no " << std::setfill(' ') << std::setw(5) << int(idx)
                  << " of approx " << N << " jobs/thread/rank." << std::endl;
      }
    }

    chi_kw = mpi::all_reduce(chi_kw);

    t.stop();
    if (c.rank() == 0) std::cout << "BSE TIME: " << double(t) << " s" << std::endl;

    return chi_kw;
  }

  // ----------------------------------------------------

  gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu_from_g_wk_and_gamma_PH(gk_iw_t g_wk, g2_iw_vt gamma_ph_wnn,
                                                                                                     int tail_corr_nwf) {

    auto _ = all_t{};

    auto target                 = gamma_ph_wnn.target();
    auto [fmesh_large, kmesh]   = g_wk.mesh();
    auto [bmesh, fmesh, fmesh2] = gamma_ph_wnn.mesh();

    double beta = fmesh.domain().beta;

    auto chi_kw = make_gf<cartesian_product<brillouin_zone, imfreq>>({kmesh, bmesh}, target);

    auto chi0_n  = make_gf<imfreq>(fmesh, target);
    auto chi0_nn = make_gf<g2_nn_t::mesh_t::var_t>({fmesh, fmesh}, target);
    auto I       = identity<Channel_t::PH>(chi0_nn);

    int nb = gamma_ph_wnn.target_shape()[0];

    gf_mesh<imfreq> fmesh_tail;
    if (tail_corr_nwf > 0)
      fmesh_tail = gf_mesh<imfreq>{beta, Fermion, tail_corr_nwf};
    else
      fmesh_tail = fmesh;

    if (fmesh_tail.size() < fmesh.size()) TRIQS_RUNTIME_ERROR << "BSE: tail size has to be larger than gamma fermi mesh.\n";

    array<std::complex<double>, 4> tr_chi(gamma_ph_wnn.target_shape());
    array<std::complex<double>, 4> tr_chi0(gamma_ph_wnn.target_shape());
    array<std::complex<double>, 4> tr_chi0_tail_corr(gamma_ph_wnn.target_shape());

    for (auto const &[k, w] : mpi_view(gf_mesh{kmesh, bmesh})) {

      triqs::utility::timer t_chi0_n, t_chi0_tr, t_bse_1, t_bse_2, t_bse_3;

      // ----------------------------------------------------
      // Build the bare bubble at k, w

      t_chi0_n.start();
      std::cout << "BSE: chi0_n ";

      auto chi0_n_tail = chi0_n_from_g_wk_PH(w, k, fmesh_tail, g_wk);
      for (auto const &n : fmesh) chi0_n[n] = chi0_n_tail(n);

      std::cout << double(t_chi0_n) << " s\n";

      // ----------------------------------------------------
      // trace the bare bubble with and without tail corrections

      t_chi0_tr.start();
      std::cout << "BSE: Tr[chi0_n] ";

      // tr_chi0_tail_corr(a, b, c, d) << density(slice_target_to_scalar(chi0_n,
      // a, b, c, d)) / beta; // does not compile

      for (auto a : range(nb)) {
        for (auto b : range(nb)) {
          for (auto c : range(nb)) {
            for (auto d : range(nb)) { tr_chi0_tail_corr(a, b, c, d) = density(slice_target_to_scalar(chi0_n_tail, a, b, c, d)) / beta; }
          }
        }
      }

      tr_chi0(a, b, c, d) << sum(chi0_n(inu)(a, b, c, d), inu = fmesh) / (beta * beta);

      std::cout << double(t_chi0_tr) << " s\n";

      // ----------------------------------------------------
      // Make two frequency object

      t_bse_1.start();
      std::cout << "BSE: chi0_nn ";

      for (auto const &n : fmesh) chi0_nn[n, n] = chi0_n[n];

      std::cout << double(t_bse_1) << " s\n";

      // ----------------------------------------------------

      t_bse_2.start();
      std::cout << "BSE: I - chi0 * gamma ";

      g2_nn_t denom = I - product<Channel_t::PH>(chi0_nn, gamma_ph_wnn[w, _, _]);

      std::cout << double(t_bse_2) << " s\n";

      t_bse_3.start();
      std::cout << "BSE: chi = [I - chi0 * gamma]^{-1} chi0 ";

      g2_nn_t chi_nn = product<Channel_t::PH>(inverse<Channel_t::PH>(denom), chi0_nn);

      std::cout << double(t_bse_3) << " s\n";

      // trace out fermionic frequencies
      std::cout << "BSE: Tr[chi] \n";
      tr_chi *= 0.0;
      for (auto const &[n1, n2] : chi_nn.mesh()) tr_chi += chi_nn[n1, n2];
      tr_chi /= beta * beta;

      // 0th order high frequency correction using the bare bubble chi0
      tr_chi += tr_chi0_tail_corr - tr_chi0;

      chi_kw[k, w] = tr_chi;
    }

    chi_kw = mpi::all_reduce(chi_kw);

    return chi_kw;
  }

  gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>>
  chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH(double mu, ek_vt e_k, g_iw_vt sigma_w, g2_iw_vt gamma_ph_wnn, int tail_corr_nwf) {

    auto _ = all_t{};

    auto target = gamma_ph_wnn.target();
    // auto [fmesh_large, kmesh] = g_wk.mesh();

    auto kmesh       = e_k.mesh();
    auto fmesh_large = sigma_w.mesh();

    auto [bmesh, fmesh, fmesh2] = gamma_ph_wnn.mesh();

    double beta = fmesh.domain().beta;

    auto chi_kw = make_gf<cartesian_product<brillouin_zone, imfreq>>({kmesh, bmesh}, target);

    auto chi0_n  = make_gf<imfreq>(fmesh, target);
    auto chi0_nn = make_gf<g2_nn_t::mesh_t::var_t>({fmesh, fmesh}, target);
    auto I       = identity<Channel_t::PH>(chi0_nn);

    int nb = gamma_ph_wnn.target_shape()[0];

    gf_mesh<imfreq> fmesh_tail;
    if (tail_corr_nwf > 0)
      fmesh_tail = gf_mesh<imfreq>{beta, Fermion, tail_corr_nwf};
    else
      fmesh_tail = fmesh;

    if (fmesh_tail.size() < fmesh.size()) TRIQS_RUNTIME_ERROR << "BSE: tail size has to be larger than gamma fermi mesh.\n";

    array<std::complex<double>, 4> tr_chi(gamma_ph_wnn.target_shape());
    array<std::complex<double>, 4> tr_chi0(gamma_ph_wnn.target_shape());
    array<std::complex<double>, 4> tr_chi0_tail_corr(gamma_ph_wnn.target_shape());

    for (auto const &[k, w] : mpi_view(gf_mesh{kmesh, bmesh})) {

      triqs::utility::timer t_chi0_n, t_chi0_tr, t_bse_1, t_bse_2, t_bse_3;

      // ----------------------------------------------------
      // Build the bare bubble at k, w

      t_chi0_n.start();
      std::cout << "BSE: chi0_n ";

      // auto chi0_n_tail = chi0_n_from_g_wk_PH(w, k, fmesh_tail, g_wk);

      auto chi0_n_tail = chi0_n_from_e_k_sigma_w_PH(w, k, fmesh_tail, mu, e_k, sigma_w);

      for (auto const &n : fmesh) chi0_n[n] = chi0_n_tail(n);

      std::cout << double(t_chi0_n) << " s\n";

      // ----------------------------------------------------
      // trace the bare bubble with and without tail corrections

      t_chi0_tr.start();
      std::cout << "BSE: Tr[chi0_n] ";

      // tr_chi0_tail_corr(a, b, c, d) << density(slice_target_to_scalar(chi0_n,
      // a, b, c, d)) / beta; // does not compile

      for (auto a : range(nb)) {
        for (auto b : range(nb)) {
          for (auto c : range(nb)) {
            for (auto d : range(nb)) { tr_chi0_tail_corr(a, b, c, d) = density(slice_target_to_scalar(chi0_n_tail, a, b, c, d)) / beta; }
          }
        }
      }

      tr_chi0(a, b, c, d) << sum(chi0_n(inu)(a, b, c, d), inu = fmesh) / (beta * beta);

      std::cout << double(t_chi0_tr) << " s\n";

      // ----------------------------------------------------
      // Make two frequency object

      t_bse_1.start();
      std::cout << "BSE: chi0_nn ";

      for (auto const &n : fmesh) chi0_nn[n, n] = chi0_n[n];

      std::cout << double(t_bse_1) << " s\n";

      // ----------------------------------------------------

      t_bse_2.start();
      std::cout << "BSE: I - chi0 * gamma ";

      g2_nn_t denom = I - product<Channel_t::PH>(chi0_nn, gamma_ph_wnn[w, _, _]);

      std::cout << double(t_bse_2) << " s\n";

      t_bse_3.start();
      std::cout << "BSE: chi = [I - chi0 * gamma]^{-1} chi0 ";

      g2_nn_t chi_nn = product<Channel_t::PH>(inverse<Channel_t::PH>(denom), chi0_nn);

      std::cout << double(t_bse_3) << " s\n";

      // trace out fermionic frequencies
      std::cout << "BSE: Tr[chi] \n";
      tr_chi *= 0.0;
      for (auto const &[n1, n2] : chi_nn.mesh()) tr_chi += chi_nn[n1, n2];
      tr_chi /= beta * beta;

      // 0th order high frequency correction using the bare bubble chi0
      tr_chi += tr_chi0_tail_corr - tr_chi0;

      chi_kw[k, w] = tr_chi;
    }

    chi_kw = mpi::all_reduce(chi_kw);

    return chi_kw;
  }

  gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu(chiq_t chiq) {

    auto mesh_k = std::get<0>(chiq.mesh());
    auto mesh_b = std::get<1>(chiq.mesh());
    auto mesh_f = std::get<2>(chiq.mesh());
    auto chiq_w = make_gf<cartesian_product<brillouin_zone, imfreq>>({mesh_k, mesh_b}, chiq.target());

    // Does not compile due to treatment of the tail (singularity)
    // chiq_w(k, iw) << sum(chiq(k, iw, inu, inup), inu=mesh, inup=mesh);

    for (auto const &[k, w, n1, n2] : chiq.mesh()) chiq_w[k, w] += chiq[k, w, n1, n2];

    double beta = mesh_f.domain().beta;
    chiq_w      = chiq_w / (beta * beta);

    return chiq_w;
  }

  gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(chiq_t chiq) {

    auto mesh_k = std::get<0>(chiq.mesh());
    auto mesh_b = std::get<1>(chiq.mesh());
    auto mesh_f = std::get<2>(chiq.mesh());
    auto chi_w  = make_gf<imfreq>(mesh_b, chiq.target());

    for (auto const &[k, w, n1, n2] : chiq.mesh()) chi_w[w] += chiq[k, w, n1, n2];

    double nk   = mesh_k.size();
    double beta = mesh_f.domain().beta;
    chi_w       = chi_w / nk / (beta * beta);

    return chi_w;
  }

} // namespace triqs_tprf
