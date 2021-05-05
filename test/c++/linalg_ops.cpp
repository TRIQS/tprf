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

#include <nda/nda.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/mc_tools/random_generator.hpp>

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace nda;

#include <triqs_tprf/linalg.hpp>

using namespace triqs_tprf;

// ----------------------------------------------------
// Helper functions

// Fill two-particle Green's function with random complex numbers Re/Im in [-1,
// 1]
void random_fill(g2_iw_t &G, int rng_seed = 23432) {
  triqs::mc_tools::random_generator RNG("mt19937", rng_seed);
  for (auto &v : G.data()) { v = std::complex<double>(RNG(2.) - 1., RNG(2.) - 1.); }
}

void print(g2_iw_t const &G) {
  typedef size_t st;
  auto data = G.data();
  foreach (data, [&data](st n1, st n2, st n3, st i, st j, st k, st l) {
    std::cout << "data(" << n1 << ", " << n2 << ", " << n3 << "; " << i << ", "
              << j << ", " << k << ", " << l
              << ") = " << data(n1, n2, n3, i, j, k, l) << "\n";
  })
    ;
}

// ----------------------------------------------------

// TESTS

// ----------------------------------------------------
// Test product with identity operator I = I * A = A * I

template <Channel_t CH> void identity_product_impl() {
  double beta = 1.2345;
  int nwb = 4, nwf = 6, norb = 3;

  gf_mesh<imfreq> fmesh{beta, Boson, nwb};
  gf_mesh<imfreq> bmesh{beta, Fermion, nwf};
  g2_iw_t A{{bmesh, fmesh, fmesh}, {norb, norb, norb, norb}};

  random_fill(A);
  // print(A);

  auto I = identity<CH>(A);

  // Test identity and product I * A = A
  {
    auto I_A = product<CH>(I, A);
    EXPECT_GF_NEAR(I_A, A);
  }

  // Test identity and product A * I = A
  {
    auto A_I = product<CH>(A, I);
    EXPECT_GF_NEAR(A_I, A);
  }
}

TEST(CtHyb, ph_identity_product) { identity_product_impl<Channel_t::PH>(); }
TEST(CtHyb, pp_identity_product) { identity_product_impl<Channel_t::PP>(); }
TEST(CtHyb, ph_bar_identity_product) {
  identity_product_impl<Channel_t::PH_bar>();
}

// ----------------------------------------------------
// Test product with inverse operator I = A^{-1} * A = A * A^{-1}

template <Channel_t CH> void inverse_product_impl() {

  double beta = 1.2345;
  int nwb = 4, nwf = 6, norb = 3;

  gf_mesh<imfreq> fmesh{beta, Boson, nwb};
  gf_mesh<imfreq> bmesh{beta, Fermion, nwf};
  g2_iw_t A{{bmesh, fmesh, fmesh}, {norb, norb, norb, norb}};

  random_fill(A);

  auto A_inv = inverse<CH>(A);
  auto I = identity<CH>(A);

  // Test A * A_inv = I
  {
    auto A_A_inv = product<CH>(A, A_inv);
    EXPECT_GF_NEAR(A_A_inv, I);
  }

  // Test identity and product A_inv * A = I
  {
    auto A_inv_A = product<CH>(A_inv, A);
    EXPECT_GF_NEAR(A_inv_A, I);
  }
}

TEST(CtHyb, ph_inverse_product) { inverse_product_impl<Channel_t::PH>(); }
TEST(CtHyb, pp_inverse_product) { inverse_product_impl<Channel_t::PP>(); }
TEST(CtHyb, ph_bar_inverse_product) {
  inverse_product_impl<Channel_t::PH_bar>();
}

// ----------------------------------------------------
MAKE_MAIN;
