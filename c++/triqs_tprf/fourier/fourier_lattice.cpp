// Copyright (c) 2014-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2014-2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2019 Simons Foundation
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
// Authors: Michel Ferrero, Olivier Parcollet, Nils Wentzell, tayral

#include "./fourier_common.hpp"
#include <itertools/itertools.hpp>

#include <fftw3.h>

namespace triqs_tprf::fourier {

  // The implementation is almost the same in both cases...
  template <typename M1, typename M2> fourier_plan __impl_plan(int fftw_backward_forward, gf_mesh<M1> const &out_mesh, gf_vec_cvt<M2> g_in) {

    //check periodization_matrix is diagonal
    auto &period_mat = g_in.mesh().periodization_matrix;
    for (auto [i, j] : itertools::product_range(period_mat.shape()[0], period_mat.shape()[1]))
      if (i != j and period_mat(i, j) != 0) {
        std::cerr
           << "WARNING: Fourier Transform of k-mesh with non-diagonal periodization matrix. Please make sure that the order of real and reciprocal space vectors is compatible for FFTW to work. (Cf. discussion doi:10.3929/ethz-a-010657714, p.26)\n";
        break;
      }

    auto g_out    = gf_vec_t<M1>{out_mesh, std::array{g_in.target_shape()[0]}};
    long n_others = g_in.data().shape()[1];

    auto dims     = g_in.mesh().get_dimensions();
    auto dims_int = stdutil::make_std_array<int>(dims);

    return _fourier_base_plan(g_in.data(), g_out.data(), dims.size(), dims_int.data(), n_others, fftw_backward_forward);
  }

  fourier_plan _fourier_plan(triqs::mesh::cyclat const &r_mesh, gf_vec_cvt<brzone> gk) { return __impl_plan(FFTW_FORWARD, r_mesh, gk); }
  fourier_plan _fourier_plan(triqs::mesh::brzone const &k_mesh, gf_vec_cvt<cyclat> gr) { return __impl_plan(FFTW_BACKWARD, k_mesh, gr); }

  // ------------------------ DIRECT TRANSFORM --------------------------------------------

  gf_vec_t<cyclat> _fourier_impl(triqs::mesh::cyclat const &r_mesh, gf_vec_cvt<brzone> gk, fourier_plan &p) {

    auto gr = gf_vec_t<cyclat>{r_mesh, {gk.target_shape()[0]}};
    _fourier_base(gk.data(), gr.data(), p);

    gr.data() /= gk.mesh().size();
    return gr;
  }

  // ------------------------ INVERSE TRANSFORM --------------------------------------------

  gf_vec_t<brzone> _fourier_impl(triqs::mesh::brzone const &k_mesh, gf_vec_cvt<cyclat> gr, fourier_plan &p) {
    auto gk = gf_vec_t<brzone>{k_mesh, {gr.target_shape()[0]}};
    _fourier_base(gr.data(), gk.data(), p);
    return std::move(gk);
  }

} // namespace triqs_tprf::fourier
