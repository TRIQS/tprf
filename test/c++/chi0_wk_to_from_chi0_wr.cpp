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

#include <triqs/clef.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs::gfs;
using namespace triqs::arrays;
using namespace triqs::lattice;

#include <triqs_tprf/types.hpp>
#include <triqs_tprf/lattice.hpp>

using namespace triqs_tprf;

TEST(lattice, chi0_wk_to_from_chi0_wr) {

 double beta = 100.0;
 int n_iw = 1028;

 int nw = 8;
 int nnu = 128;

 int nk = 4; 
 double t = 1.0;
 auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
 
 triqs::clef::placeholder<0> om_;
 triqs::clef::placeholder<1> k_;

 auto e_k = ek_t{{bz, nk}, {1, 1}};
 e_k(k_) << - 2*t * (cos(k_(0)) + cos(k_(1)));

 double mu = 0.;
 auto mesh = g_iw_t::mesh_t{beta, Fermion, n_iw};

 auto sigma_w = g_iw_t{mesh, {1, 1}};
 sigma_w(om_) << 1./om_;
 
 auto g_wk = lattice_dyson_g_wk(mu, e_k, sigma_w);
 auto g_wr = fourier_wk_to_wr(g_wk);
 
 auto chi0_wr = chi0r_from_gr_PH(nw, nnu, g_wr);
 auto chi0_wk = chi0q_from_chi0r(chi0_wr);
 auto chi0_wr_ref = chi0r_from_chi0q(chi0_wk);

 EXPECT_ARRAY_NEAR(chi0_wr.data(), chi0_wr_ref.data()); 
}

MAKE_MAIN;
