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

#include "types.hpp"
#include "lattice.hpp"
using namespace tprf;

TEST(lattice, g0k_to_from_g0r) {
 double beta = 100.0;
 int n_iw = 1025;

 int nk = 4; 
 double t = 1.0;
 auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
 
 triqs::clef::placeholder<0> om_;
 triqs::clef::placeholder<1> k_;

 auto ek = ek_t{{bz, nk}, {1, 1}};
 ek(k_) << - 2*t * (cos(k_(0)) + cos(k_(1)));

 double mu = 0.;
 auto mesh = g_iw_t::mesh_t{beta, Fermion, n_iw};
 auto g0k = g0k_from_ek(mu, ek, mesh);

 auto g0r = gr_from_gk(g0k);
 auto g0k_ref = gk_from_gr(g0r);

 EXPECT_CLOSE_ARRAY(g0k.data(), g0k_ref.data()); 
}

TEST(lattice, gk_to_from_gr) {
 double beta = 100.0;
 int n_iw = 1025;

 int nk = 4; 
 double t = 1.0;
 auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
 
 triqs::clef::placeholder<0> om_;
 triqs::clef::placeholder<1> k_;

 auto ek = ek_t{{bz, nk}, {1, 1}};
 ek(k_) << - 2*t * (cos(k_(0)) + cos(k_(1)));

 double mu = 0.;
 auto mesh = g_iw_t::mesh_t{beta, Fermion, n_iw};

 auto sigma = g_iw_t{mesh, {1, 1}};
 sigma(om_) << 1./om_;
 
 auto gk = gk_from_ek_sigma(mu, ek, sigma);

 auto gr = gr_from_gk(gk);
 auto gk_ref = gk_from_gr(gr);

 EXPECT_CLOSE_ARRAY(gk.data(), gk_ref.data()); 
}

MAKE_MAIN;
