// Copyright (c) 2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2020 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Hugo Strand, Michel Ferrero, Nils Wentzell

#pragma once

#include <nda/nda.hpp>

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/utility/tuple_tools.hpp>

#include "fourier_common.hpp"

namespace triqs_tprf::fourier {

  using namespace nda;
  using namespace triqs::gfs;
  using namespace triqs::mesh;

  /*------------------------------------------------------------------------------------------------------
            Implementation
  *-----------------------------------------------------------------------------------------------------*/

  template <typename V> using gf_vec_t   = gf<V, tensor_valued<1>>;
  template <typename V> using gf_vec_vt  = gf_view<V, tensor_valued<1>>;
  template <typename V> using gf_vec_cvt = gf_const_view<V, tensor_valued<1>>;

  // matsubara
  gf_vec_t<imfreq> _fourier_impl(mesh::imfreq const &iw_mesh, gf_vec_cvt<imtime> gt, fourier_plan &p, array_const_view<dcomplex, 2> mom_23 = {});
  gf_vec_t<imtime> _fourier_impl(mesh::imtime const &tau_mesh, gf_vec_cvt<imfreq> gw, fourier_plan &p, array_const_view<dcomplex, 2> mom_123 = {});
  fourier_plan _fourier_plan(mesh::imfreq const &iw_mesh, gf_vec_cvt<imtime> gt);
  fourier_plan _fourier_plan(mesh::imtime const &tau_mesh, gf_vec_cvt<imfreq> gw);

  // lattice
  gf_vec_t<cyclat> _fourier_impl(mesh::cyclat const &r_mesh, gf_vec_cvt<brzone> gk, fourier_plan &p);
  gf_vec_t<brzone> _fourier_impl(mesh::brzone const &k_mesh, gf_vec_cvt<cyclat> gr, fourier_plan &p);
  fourier_plan _fourier_plan(mesh::cyclat const &r_mesh, gf_vec_cvt<brzone> gk);
  fourier_plan _fourier_plan(mesh::brzone const &k_mesh, gf_vec_cvt<cyclat> gr);

  /*------------------------------------------------------------------------------------------------------
   *
   * The general Fourier function
   * gin : input
   * gout : output
   * opt_args : e.g. moments. Must be flatten_2d
   *
   *-----------------------------------------------------------------------------------------------------*/

  // this function just regroups the function, and calls the vector_valued gf core
  // implementation
  template <int N, typename V1, typename V2, typename T, typename... OptArgs>
  void _fourier_with_plan(gf_const_view<V1, T> gin, gf_view<V2, T> gout, fourier_plan &p, OptArgs const &...opt_args) {

    // pb std::get<0> would not work on a non composite mesh. We use a little lambda to deduce ref and type
    auto const &out_mesh = [&gout]() -> auto const & { // NB must return a reference
      using m_t = std::decay_t<decltype(gout.mesh())>;
      if constexpr (triqs::mesh::is_product<m_t>)
        return std::get<N>(gout.mesh());
      else
        return gout.mesh();
    }
    ();

    auto gout_flatten = _fourier_impl(out_mesh, flatten_gf_2d<N>(gin), p, flatten_2d<0>(make_array_const_view(opt_args))...);
    auto _            = ellipsis();
    if constexpr (gin.data_rank == 1)
      gout.data() = gout_flatten.data()(_, 0); // gout is scalar, gout_flatten vectorial
    else {
      // inverse operation as flatten_2d, exactly
      auto g_rot = nda::rotate_index_view<N>(gout.data());
      for (auto mp : out_mesh) {
        auto g_rot_sl = g_rot(mp.data_index(), _); // if the array is long, it is faster to precompute the view ...
        auto gout_col = gout_flatten.data()(mp.data_index(), _);
        nda::for_each(g_rot_sl.shape(), [&g_rot_sl, &gout_col, c = long(0)](auto &&...i) mutable { return g_rot_sl(i...) = gout_col(c++); });
      }
    }
  }

  // this function just regroups the function, and calls the vector_valued gf core
  // implementation
  template <int N, typename V1, typename V2, typename T, typename... OptArgs>
  fourier_plan _fourier_plan(gf_const_view<V1, T> gin, gf_view<V2, T> gout, OptArgs const &...opt_args) {

    // pb std::get<0> would not work on a non composite mesh. We use a little lambda to deduce ref and type
    auto const &out_mesh = [&gout]() -> auto const & { // NB must return a reference
      using m_t = std::decay_t<decltype(gout.mesh())>;
      if constexpr (triqs::mesh::is_product<m_t>)
        return std::get<N>(gout.mesh());
      else
        return gout.mesh();
    }
    ();

    return _fourier_plan(out_mesh, flatten_gf_2d<N>(gin), flatten_2d(opt_args, N)...);
  }
} // namespace triqs_tprf::fourier
