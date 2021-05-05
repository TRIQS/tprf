//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include <triqs/test_tools/gfs.hpp>

using namespace nda::clef;
using namespace triqs::lattice;

int n_k     = 4; // 16; // too long, 30 s
int nw      = 3;
double beta = 20, mu = 0.;

placeholder<0> k_;
placeholder<1> q_;
placeholder<2> r_;
placeholder<3> iw_;
placeholder<4> inu_;

placeholder<5> a;
placeholder<6> b;
placeholder<7> c;
placeholder<8> d;


auto _ = all_t{};

auto bz     = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
auto eps_k_ = -2 * (cos(k_(0)) + cos(k_(1)));

// ------------------------------------------------------------

TEST(Gf, Bubble) {
 auto chi0q = gf<prod<imfreq, imfreq, brillouin_zone>>{{{beta, Fermion, nw}, {beta, Boson, nw}, {bz, n_k}}, {1, 1}};
 auto chi0r = gf<prod<imfreq, imfreq, cyclat>>{{{beta, Fermion, nw}, {beta, Boson, nw}, {n_k, n_k}}, {1, 1}};

 auto chi0q_from_r = chi0q;

 auto Gk = gf<prod<imfreq, brillouin_zone>, matrix_valued>{{{beta, Fermion, nw}, {bz, n_k}}, {1, 1}};
 auto Gr = gf<prod<imfreq, cyclat>, matrix_valued>{{{beta, Fermion, nw}, {n_k, n_k}}, {1, 1}};

 Gk(inu_, k_) << 1 / (inu_ + mu - eps_k_);

 auto k_mesh = std::get<1>(Gk.mesh());

 chi0q(inu_, iw_, q_) << sum(Gk(inu_, k_) * Gk(inu_ + iw_, k_ + q_), k_ = k_mesh) / k_mesh.size();

 for (auto const &inu : std::get<0>(Gr.mesh()))
  Gr[inu, _] = fourier(Gk[inu, _]);

 chi0r(inu_, iw_, r_) << Gr(inu_, r_) * Gr(inu_ + iw_, -r_);

 for (auto const &inu : std::get<0>(chi0q_from_r.mesh())) {
  for (auto const &iw : std::get<1>(chi0q_from_r.mesh())) {

   auto size_from_r = chi0q_from_r[inu, iw, _].mesh().size();
   auto size_r = chi0r[inu, iw, _].mesh().size();
     
   EXPECT_EQ(size_from_r, 16);
   EXPECT_EQ(size_r, 16);

   chi0q_from_r[inu, iw, _] = fourier(chi0r[inu, iw, _]);
  }
 }
 EXPECT_ARRAY_NEAR(chi0q_from_r.data(), chi0q.data());
}

// ------------------------------------------------------------

TEST(Gf, BubbleScalar) {

 auto chi0q =
     gf<prod<imfreq, imfreq, brillouin_zone>, scalar_valued>{{{beta, Fermion, nw}, {beta, Boson, nw}, {bz, n_k}}};
 auto chi0r =
     gf<prod<imfreq, imfreq, cyclat>, scalar_valued>{{{beta, Fermion, nw}, {beta, Boson, nw}, {n_k, n_k}}};

 auto chi0q_from_r = chi0q;

 auto Gk = gf<prod<imfreq, brillouin_zone>, scalar_valued>{{{beta, Fermion, nw}, {bz, n_k}}};
 auto Gr = gf<prod<imfreq, cyclat>, scalar_valued>{{{beta, Fermion, nw}, {n_k, n_k}}};

 Gk(inu_, k_) << 1 / (inu_ + mu - eps_k_);

 auto Gmesh = std::get<1>(Gk.mesh());
 chi0q(inu_, iw_, q_) << sum(Gk(inu_, k_) * Gk(inu_ + iw_, k_ + q_), k_ = Gmesh) / Gmesh.size();

 for (auto const &inu : std::get<0>(Gr.mesh()))
  Gr[inu, _] = fourier(Gk[inu, _]);

 chi0r(inu_, iw_, r_) << Gr(inu_, r_) * Gr(inu_ + iw_, -r_);

 for (auto const &inu : std::get<0>(chi0q_from_r.mesh())) {
  for (auto const &iw : std::get<1>(chi0q_from_r.mesh())) {
   chi0q_from_r[inu, iw, _] = fourier(chi0r[inu, iw, _]);
  }
 }
 EXPECT_ARRAY_NEAR(chi0q_from_r.data(), chi0q.data());
}

// ------------------------------------------------------------

TEST(Gf, BubbleTensor) {

 auto chi0q =
   gf<prod<imfreq, imfreq, brillouin_zone>, tensor_valued<4>>{{{beta, Fermion, nw}, {beta, Boson, nw}, {bz, n_k}}, {1, 1, 1, 1}};
 auto chi0r =
   gf<prod<imfreq, imfreq, cyclat>, tensor_valued<4>>{{{beta, Fermion, nw}, {beta, Boson, nw}, {n_k, n_k}}, {1, 1, 1, 1}};

 auto chi0q_from_r = chi0q;

 auto Gk = gf<prod<imfreq, brillouin_zone>, matrix_valued>{{{beta, Fermion, nw}, {bz, n_k}}, {1, 1}};
 auto Gr = gf<prod<imfreq, cyclat>, matrix_valued>{{{beta, Fermion, nw}, {n_k, n_k}}, {1, 1}};

 Gk(inu_, k_)(a, b) << kronecker(a, b) * 1 / (inu_ + mu - eps_k_);

 auto Gmesh = std::get<1>(Gk.mesh());

 // PH chi0q
 chi0q(inu_, iw_, q_)(a, b, c, d) << sum(Gk(inu_, k_)(d, a) * Gk(inu_ + iw_, k_ + q_)(b, c), k_ = Gmesh) / Gmesh.size();

 for (auto const &inu : std::get<0>(Gr.mesh()))
  Gr[inu, _] = fourier(Gk[inu, _]);

 chi0r(inu_, iw_, r_)(a, b, c, d) << Gr(inu_, r_)(d, a) * Gr(inu_ + iw_, -r_)(b, c);

 for (auto const &inu : std::get<0>(chi0q_from_r.mesh())) {
  for (auto const &iw : std::get<1>(chi0q_from_r.mesh())) {
   chi0q_from_r[inu, iw, _] = fourier(chi0r[inu, iw, _]);
  }
 }
 EXPECT_ARRAY_NEAR(chi0q_from_r.data(), chi0q.data());
}

// ------------------------------------------------------------

MAKE_MAIN;
