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

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace nda;
using namespace triqs::lattice;

TEST(tprf, gf_inverse) {
 double beta = 100.0;
 int n_iw = 1025;

 int nk = 4; 
 double t = 1.0;
 auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
 
 auto G_iw = gf<prod<imfreq, brzone>>{{{beta, Fermion, n_iw}, {bz, nk}}, {1, 1}};

 nda::clef::placeholder<0> om_;
 nda::clef::placeholder<1> k_;
 
 G_iw(om_, k_) << om_ - 2*t * (cos(k_(0)) + cos(k_(1)));

 //G_iw = inverse(G_iw); // does not work, see TRIQS issue #463

 for (auto kidx : std::get<1>(G_iw.mesh())) {
   auto _ = all_t{};
   auto G = G_iw[_, kidx];
   G_iw[_, kidx] = inverse(G);
 }
}

MAKE_MAIN;
