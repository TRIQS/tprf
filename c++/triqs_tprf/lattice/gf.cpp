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

#include <nda/linalg.hpp>
using nda::inverse;

#include "gf.hpp"

#include <omp.h>
#include "../mpi.hpp"
#include "fourier.hpp"

namespace triqs_tprf {

  namespace {
    using namespace fourier;
  }
 

// ----------------------------------------------------
// g

g_wk_t lattice_dyson_g0_wk(double mu, e_k_cvt e_k, gf_mesh<imfreq> mesh, mpi::communicator const &c) {

  auto I = nda::eye<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_wk_t g0_wk({mesh, e_k.mesh()}, e_k.target_shape());
  g0_wk() = 0.0;

  auto arr = mpi_view(g0_wk.mesh(), c);

#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);
    g0_wk[w, k] = inverse((w + mu)*I - e_k(k));      
  }

  g0_wk = mpi::all_reduce(g0_wk, c);
  return g0_wk;
}

// ----------------------------------------------------
// g0 real frequencies

g_fk_t lattice_dyson_g0_fk(double mu, e_k_cvt e_k, gf_mesh<refreq> mesh, double delta) {

  auto I = nda::eye<ek_vt::scalar_t>(e_k.target_shape()[0]);
  g_fk_t g0_fk({mesh, e_k.mesh()}, e_k.target_shape());
  g0_fk() = 0.0;
  std::complex<double> idelta(0.0, delta);

  auto arr = mpi_view(g0_fk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[f, k] = arr(idx);
    g0_fk[f, k]  = inverse((f + idelta + mu) * I - e_k(k));
  }

  g0_fk = mpi::all_reduce(g0_fk);
  return g0_fk;
}

// ----------------------------------------------------

template<typename sigma_t>
auto lattice_dyson_g_generic(double mu, e_k_cvt e_k, sigma_t sigma, mpi::communicator const &c){

  auto &freqmesh = [&sigma]() -> auto & {
    if constexpr (sigma_t::arity == 1) return sigma.mesh();
    else return std::get<0>(sigma.mesh());
  }();

  using scalar_t = e_k_cvt::scalar_t;
  auto I = nda::eye<scalar_t>(e_k.target_shape()[0]);
  
  g_wk_t g_wk({freqmesh, e_k.mesh()}, e_k.target_shape());
  g_wk() = 0.0;

  auto arr = mpi_view(g_wk.mesh(), c);
#pragma omp parallel for
  for (unsigned int idx = 0; idx < arr.size(); idx++) {
    auto &[w, k] = arr(idx);

    array<scalar_t, 2> sigmaterm;
    if constexpr (sigma_t::arity == 1) sigmaterm = sigma[w];
    else sigmaterm = sigma[w, k];

    g_wk[w, k] = inverse((w + mu)*I - e_k(k) - sigmaterm);
  }

  g_wk = mpi::all_reduce(g_wk, c);
  return g_wk;
}


g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_wk_cvt sigma_wk, mpi::communicator const &c) {
  return lattice_dyson_g_generic(mu, e_k, sigma_wk, c);
}
  

g_wk_t lattice_dyson_g_wk(double mu, e_k_cvt e_k, g_w_cvt sigma_w, mpi::communicator const &c) {
  return lattice_dyson_g_generic(mu, e_k, sigma_w, c);
}


// ----------------------------------------------------
// g in real frequencies

g_fk_t lattice_dyson_g_fk(double mu, e_k_cvt e_k, g_fk_cvt sigma_fk, double delta) {

  auto I    = nda::eye<ek_vt::scalar_t>(e_k.target_shape()[0]);
  auto g_fk = make_gf(sigma_fk);
  g_fk()    = 0.0;
  std::complex<double> idelta(0.0, delta);

  auto arr = mpi_view(g_fk.mesh());
#pragma omp parallel for
  for (int idx = 0; idx < arr.size(); idx++) {
    auto &[f, k] = arr(idx);
    g_fk[f, k]   = inverse((f + idelta + mu) * I - e_k(k) - sigma_fk[f, k]);
  }

  g_fk = mpi::all_reduce(g_fk);

  return g_fk;
}

// ----------------------------------------------------

g_w_t lattice_dyson_g_w(double mu, e_k_cvt e_k, g_w_cvt sigma_w, mpi::communicator const &c) {

  auto g_wk = lattice_dyson_g_generic(mu, e_k, sigma_w, c);
  auto &[wmesh, kmesh] = g_wk.mesh();

  g_w_t g_w(wmesh, e_k.target_shape());
  g_w() = 0.0;

  for (auto const &[w, k] : mpi_view(g_wk.mesh(), c))
    g_w[w] += g_wk[w, k];

  g_w = mpi::all_reduce(g_w, c);
  g_w /= kmesh.size();
  return g_w;
}

// ----------------------------------------------------
// Transformations: real space <-> reciprocal space 
  
g_wr_t fourier_wk_to_wr(g_wk_cvt g_wk, mpi::communicator const &c) {
  auto g_wr = fourier_wk_to_wr_general_target(g_wk, c);
  return g_wr;
}

g_wk_t fourier_wr_to_wk(g_wr_cvt g_wr, mpi::communicator const &c) {
  auto g_wk = fourier_wr_to_wk_general_target(g_wr, c);
  return g_wk;
}

g_fr_t fourier_fk_to_fr(g_fk_cvt g_fk) {
  auto g_fr = fourier_wk_to_wr_general_target(g_fk);
  return g_fr;
}

g_fk_t fourier_fr_to_fk(g_fr_cvt g_fr) {
  auto g_fk = fourier_wr_to_wk_general_target(g_fr);
  return g_fk;
}

g_Tr_t fourier_Tk_to_Tr(g_Tk_cvt g_Tk) {
  auto g_Tr = fourier_wk_to_wr_general_target(g_Tk);
  return g_Tr;
}

g_Tk_t fourier_Tr_to_Tk(g_Tr_cvt g_Tr) {
  auto g_Tk = fourier_wr_to_wk_general_target(g_Tr);
  return g_Tk;
}

chi_Tk_t fourier_Tr_to_Tk(chi_Tr_cvt chi_Tr) {
  auto chi_Tk = fourier_wr_to_wk_general_target(chi_Tr);
  return chi_Tk;
}

// ----------------------------------------------------
// Transformations: Matsubara frequency <-> imaginary time

g_wr_t fourier_tr_to_wr(g_tr_cvt g_tr, int nw, mpi::communicator const &c) {
  auto g_wr = fourier_tr_to_wr_general_target(g_tr, nw, c);
  return g_wr;
}

g_tr_t fourier_wr_to_tr(g_wr_cvt g_wr, int nt, mpi::communicator const &c) {
  auto g_tr = fourier_wr_to_tr_general_target(g_wr, nt, c);
  return g_tr;
}

} // namespace triqs_tprf
