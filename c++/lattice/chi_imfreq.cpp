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

#include "../linalg.hpp"

#include "common.hpp"
#include "chi_imfreq.hpp"

namespace tprf {

// ----------------------------------------------------
// chi0 bubble in Matsubara frequency

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

chi0r_t chi0r_from_chi0q(chi0q_vt chi0q) {

  auto mb = std::get<0>(chi0q.mesh());
  auto mf = std::get<1>(chi0q.mesh());
  auto bzmesh = std::get<2>(chi0q.mesh());

  auto clmesh = gf_mesh<cyclic_lattice>{bzmesh.domain().lattice(), bzmesh.periodization_matrix};
  
  auto chi0r =
      make_gf<chi0r_t::mesh_t::var_t>({mb, mf, clmesh}, chi0q.target());

  for (auto const &widx : std::get<0>(chi0r.mesh())) {
    for (auto const &nidx : std::get<1>(chi0r.mesh())) {
      auto _ = all_t{};
      //chi0r[widx, nidx, _] = inverse_fourier(chi0q[widx, nidx, _]);
      chi0r[widx, nidx, _] = fourier(chi0q[widx, nidx, _]);
    }
  }
  return chi0r;
}

chi0q_t chi0q_from_chi0r(chi0r_vt chi0r) {

  auto mb = std::get<0>(chi0r.mesh());
  auto mf = std::get<1>(chi0r.mesh());

  auto clmesh = std::get<2>(chi0r.mesh());
  auto bzmesh = gf_mesh<brillouin_zone>{brillouin_zone{clmesh.domain()}, clmesh.periodization_matrix};
  
  auto chi0q =
      make_gf<chi0q_t::mesh_t::var_t>({mb, mf, bzmesh}, chi0r.target());

  for (auto const &widx : std::get<0>(chi0q.mesh())) {
    for (auto const &nidx : std::get<1>(chi0q.mesh())) {
      auto _ = all_t{};
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

  chi0q_w(iw, k) << sum(chi0q(iw, inu, k), inu = mesh) / (beta * beta);

  return chi0q_w;
}

gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>>
chi0q_sum_nu_tail_corr_PH(chi0q_t chi0q) {

  auto mesh = std::get<1>(chi0q.mesh());
  auto chi0q_w = make_gf<cartesian_product<imfreq, brillouin_zone>>(
      {std::get<0>(chi0q.mesh()), std::get<2>(chi0q.mesh())}, chi0q.target());

  int nb = chi0q.target_shape()[0];
  double beta = mesh.domain().beta;

  auto wmesh = std::get<0>(chi0q.mesh());
  auto qmesh = std::get<2>(chi0q.mesh());

  for (auto const &w : wmesh) {
    for (auto const &q : qmesh) {
      auto _ = all_t{};

      auto chi = chi0q[w, _, q];

      for (auto a : range(nb)) {
        for (auto b : range(nb)) {
          for (auto c : range(nb)) {
            for (auto d : range(nb)) {
              auto chi_abcd = slice_target_to_scalar(chi, a, b, c, d);
              chi0q_w[w, q](a, b, c, d) = density(chi_abcd) / beta;
            }
          }
        }
      }
    }
  }
  return chi0q_w;
}

gf<imfreq, tensor_valued<4>> chi0q_sum_nu_q(chi0q_t chi0q) {

  auto mesh_b = std::get<0>(chi0q.mesh());
  auto mesh_f = std::get<1>(chi0q.mesh());
  auto mesh_k = std::get<2>(chi0q.mesh());

  auto chi0_w = make_gf<imfreq>(mesh_b, chi0q.target());

  for (auto const & [ w, n, k ] : chi0q.mesh())
    chi0_w[w] += chi0q[w, n, k];

  double nk = mesh_k.size();
  double beta = mesh_f.domain().beta;
  chi0_w = chi0_w / nk / (beta * beta);

  return chi0_w;
}

// ----------------------------------------------------
// chi

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

gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>>
chiq_sum_nu(chiq_t chiq) {

  auto mesh_k = std::get<0>(chiq.mesh());
  auto mesh_b = std::get<1>(chiq.mesh());
  auto mesh_f = std::get<2>(chiq.mesh());
  auto chiq_w = make_gf<cartesian_product<brillouin_zone, imfreq>>(
      {mesh_k, mesh_b}, chiq.target());

  // Does not compile due to treatment of the tail (singularity)
  // chiq_w(k, iw) << sum(chiq(k, iw, inu, inup), inu=mesh, inup=mesh);

  for (auto const & [ k, w, n1, n2 ] : chiq.mesh())
    chiq_w[k, w] += chiq[k, w, n1, n2];

  double beta = mesh_f.domain().beta;
  chiq_w = chiq_w / (beta * beta);

  return chiq_w;
}

gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(chiq_t chiq) {

  auto mesh_k = std::get<0>(chiq.mesh());
  auto mesh_b = std::get<1>(chiq.mesh());
  auto mesh_f = std::get<2>(chiq.mesh());
  auto chi_w = make_gf<imfreq>(mesh_b, chiq.target());

  for (auto const & [ k, w, n1, n2 ] : chiq.mesh())
    chi_w[w] += chiq[k, w, n1, n2];

  double nk = mesh_k.size();
  double beta = mesh_f.domain().beta;
  chi_w = chi_w / nk / (beta * beta);

  return chi_w;
}

} // namespace tprf
