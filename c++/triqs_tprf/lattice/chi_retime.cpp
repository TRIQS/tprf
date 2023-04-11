/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2022, H. U.R. Strand
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

#include <nda/linalg.hpp>

#include "common.hpp"
#include "../mpi.hpp"
#include "chi_retime.hpp"

#include "../fourier/fourier.hpp"
#include "fourier.hpp"

namespace triqs_tprf {

  namespace {
    using namespace fourier;
  }

// ----------------------------------------------------
// helper functions
  
double fermi_scalar(double e) {
  if( e < 0 ) {
    return 1. / (exp(e) + 1.);
  } else {
    double exp_me = exp(-e);
    return exp_me / (1 + exp_me);
  }
}

auto fermi_nda = nda::map([](double e) { return fermi_scalar(e); });
  
// ----------------------------------------------------
// g0_Tk_les, g0_Tk_gtr in real time

std::tuple<g_Tk_t, g_Tk_t> g0_Tk_les_gtr_from_e_k(e_k_cvt e_k, mesh::retime Tmesh, double beta) {

  auto I = std::complex(0.,1.);
  auto kmesh = e_k.mesh();
  
  g_Tk_t g0_Tk_les({Tmesh, kmesh}, e_k.target_shape());
  g_Tk_t g0_Tk_gtr({Tmesh, kmesh}, e_k.target_shape());

  auto arr = mpi_view(kmesh);

#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto & k = arr(idx);

    matrix<std::complex<double>> e_k_mat(e_k[k]);
    auto [ek, Uk] = linalg::eigenelements(e_k_mat);

    auto occ = fermi_nda(beta*ek);
    
    for (auto const &T : Tmesh) {
      auto exp_T = exp(-I * ek * double(T));
      auto occ_exp_les = +I * occ * exp_T;
      auto occ_exp_gtr = -I * (1. - occ) * exp_T;

      for (auto [a, b] : g0_Tk_les.target_indices())
        for (auto c : range(ek.size())) {
          g0_Tk_les[T, k](a, b) += Uk(a, c) * occ_exp_les(c) * dagger(Uk)(c, b);
          g0_Tk_gtr[T, k](a, b) += Uk(a, c) * occ_exp_gtr(c) * dagger(Uk)(c, b);
        }
    }
  }

  g0_Tk_les = mpi::all_reduce(g0_Tk_les);
  g0_Tk_gtr = mpi::all_reduce(g0_Tk_gtr);
  
  return {g0_Tk_les, g0_Tk_gtr};
}
  
// ----------------------------------------------------
// chi0 bubble in real time
  
chi_Tr_t chi0_Tr_from_g_Tr_PH(g_Tr_cvt g_Tr_les, g_Tr_cvt g_Tr_gtr) {

  auto _ = all_t{};
  auto I = std::complex(0.,1.);
  
  auto Tmesh = std::get<0>(g_Tr_les.mesh());
  auto rmesh = std::get<1>(g_Tr_les.mesh());

  int nb = g_Tr_les.target().shape()[0];

  chi_Tr_t chi0_Tr{{Tmesh, rmesh}, {nb, nb, nb, nb}};

  auto g_target = g_Tr_les.target();
  auto chi_target = chi0_Tr.target();
  
  //for (auto const &r : rmesh) {

  auto arr = mpi_view(Tmesh);

#pragma omp parallel for 
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto & T = arr(idx);
    
    for (auto const &r : rmesh)
      for (auto [a, b, c, d] : chi0_Tr.target_indices())
        chi0_Tr[T, r](a, b, c, d) = +I * g_Tr_les[T, r](d, a) * conj(g_Tr_gtr[T, -r](b, c)) - I * g_Tr_gtr[T, r](d, a) * conj(g_Tr_les[T, -r](b, c));
  }

  chi0_Tr = mpi::all_reduce(chi0_Tr);

  return chi0_Tr;
}
  
} // namespace triqs_tprf
  
