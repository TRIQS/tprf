
#include <triqs/test_tools/gfs.hpp>

using namespace triqs::clef;
using namespace triqs::lattice;

#include "mpi.hpp"

using namespace tprf;

// ------------------------------------------------------------

TEST(Gf, Bubble) {

  triqs::mpi::communicator c;
  
  double beta = 20;
  int n_k = 2;
  int nw = 3;

  auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
  auto g_wk = gf<cartesian_product<imfreq, brillouin_zone>>{
      {{beta, Fermion, nw}, {bz, n_k}}, {1, 1}};

  auto [wmesh, kmesh] = g_wk.mesh();
  
  // loop over frequency and momentum in mpi

  c.barrier();

  for (auto const &[w, k] : mpi_view(g_wk.mesh())) {
    std::cout << "rank " << c.rank() << " size " << c.size() << " : " << w.linear_index()
              << ", " << k.linear_index() << std::endl;
  }

  // loop over frequency in mpi

  c.barrier();

  for (auto const &w : mpi_view(wmesh)) {
    std::cout << "rank " << c.rank() << " size " << c.size() << " : " << w.linear_index() << std::endl;
  }

  // loop over momentum in mpi

  c.barrier();

  for (auto const &k : mpi_view(kmesh)) {
    std::cout << "rank " << c.rank() << " size " << c.size() << " : " << k.linear_index() << std::endl;
  }
  
  c.barrier();
}

// ------------------------------------------------------------

MAKE_MAIN;
