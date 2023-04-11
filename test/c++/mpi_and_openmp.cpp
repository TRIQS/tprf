
#include <omp.h>

#include <triqs/test_tools/gfs.hpp>

using namespace nda::clef;
using namespace triqs::lattice;

#include <triqs_tprf/mpi.hpp>

using namespace triqs_tprf;

// ------------------------------------------------------------

TEST(mpi, mpi_view) {

  mpi::communicator c;

  double beta = 20;
  int n_k = 2;
  int nw = 3;

  auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
  auto g_wk = gf<prod<imfreq, brzone>>{
      {{beta, Fermion, nw}, {bz, n_k}}, {1, 1}};

  auto [wmesh, kmesh] = g_wk.mesh();

  // loop over frequency and momentum in mpi

  c.barrier();

  for (auto [w, k] : mpi_view(g_wk.mesh())) {
    std::cout << "rank " << c.rank() << " size " << c.size() << " : " << w.data_index() << ", " << k.data_index() << std::endl;
  }

  // loop over frequency in mpi

  c.barrier();

  for (auto w : mpi_view(wmesh)) { std::cout << "rank " << c.rank() << " size " << c.size() << " : " << w.data_index() << std::endl; }

  // loop over momentum in mpi

  c.barrier();

  for (auto k : mpi_view(kmesh)) { std::cout << "rank " << c.rank() << " size " << c.size() << " : " << k.data_index() << std::endl; }

  c.barrier();
}

// ------------------------------------------------------------

TEST(mpi, mpi_view_openmp) {

  mpi::communicator c;

  double beta = 20;
  int n_k = 2;
  int nw = 3;

  auto bz = brillouin_zone{bravais_lattice{{{1, 0}, {0, 1}}}};
  auto g_wk = gf<prod<imfreq, brzone>>{
      {{beta, Fermion, nw}, {bz, n_k}}, {1, 1}};

  auto [wmesh, kmesh] = g_wk.mesh();

  // loop over frequency and momentum in mpi + openmp

  c.barrier();

  {
    auto arr = mpi_view(g_wk.mesh());
#pragma omp parallel for
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &[w, k] = arr[idx];

      int tid = omp_get_thread_num();
      std::cout << "thread " << tid << " rank " << c.rank() << " size " << c.size() << " : " << w.data_index() << ", " << k.data_index() << std::endl;
    }
  }

  // loop over frequency in mpi + openmp

  c.barrier();

  {
    auto arr = mpi_view(wmesh);
#pragma omp parallel for
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &w = arr[idx];

      int tid = omp_get_thread_num();
      std::cout << "thread " << tid << " rank " << c.rank() << " size " << c.size() << " : " << w.data_index() << std::endl;
    }
  }

  // loop over momentum in mpi + openmp

  c.barrier();

  {
    auto arr = mpi_view(kmesh);
#pragma omp parallel for
    for (unsigned int idx = 0; idx < arr.size(); idx++) {
      auto &k = arr[idx];

      int tid = omp_get_thread_num();
      std::cout << "thread " << tid << " rank " << c.rank() << " size " << c.size() << " : " << k.data_index() << std::endl;
    }
  }

  c.barrier();
}

MAKE_MAIN;
