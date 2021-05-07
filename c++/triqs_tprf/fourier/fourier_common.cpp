// Copyright (c) 2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018 Simons Foundation
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

#include "./fourier_common.hpp"

#include <fftw3.h>

namespace triqs_tprf::fourier {

  fourier_plan _fourier_base_plan(array_const_view<dcomplex, 2> in, array_view<dcomplex, 2> out, int rank, int *dims, int fftw_count,
                                  int fftw_backward_forward) {

    auto in_fft  = reinterpret_cast<fftw_complex *>(const_cast<dcomplex *>(in.data()));
    auto out_fft = reinterpret_cast<fftw_complex *>(out.data());

    auto p = fftw_plan_many_dft(rank,                        // rank
                                dims,                        // the dimension
                                fftw_count,                  // how many FFT : here 1
                                in_fft,                      // in data
                                NULL,                        // embed : unused. Doc unclear ?
                                in.indexmap().strides()[0],  // stride of the in data
                                1,                           // in : shift for multi fft.
                                out_fft,                     // out data
                                NULL,                        // embed : unused. Doc unclear ?
                                out.indexmap().strides()[0], // stride of the out data
                                1,                           // out : shift for multi fft.
                                fftw_backward_forward, FFTW_ESTIMATE);

    return {(void *)p, [](void *p) { fftw_destroy_plan((fftw_plan)p); }};
  }

  void _fourier_base(array_view<dcomplex, 2> in, array_view<dcomplex, 2> out, fourier_plan &plan) {

    auto in_fft  = reinterpret_cast<fftw_complex *>(in.data());
    auto out_fft = reinterpret_cast<fftw_complex *>(out.data());

    fftw_execute_dft((fftw_plan)plan.get(), in_fft, out_fft);
  }

} // namespace triqs_tprf::fourier
