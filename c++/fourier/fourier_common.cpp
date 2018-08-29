/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 * Copyright (C) 2018- by Simons Foundation
 *               authors : O. Parcollet, N. Wentzell, H. U.R. Strand
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
#include "./fourier_common.hpp"
#include <triqs/gfs.hpp>

namespace tprf::fourier {

void _fourier_base(array_const_view<dcomplex, 2> in,
                   array_view<dcomplex, 2> out, int rank, int *dims,
                   int fftw_count, int fftw_backward_forward) {

  auto in_fft = reinterpret_cast<fftw_complex *>(in.data_start());
  auto out_fft = reinterpret_cast<fftw_complex *>(out.data_start());

  auto p =
      fftw_plan_many_dft(rank,       // rank
                         dims,       // the dimension
                         fftw_count, // how many FFT : here 1
                         in_fft,     // in data
                         NULL,       // embed : unused. Doc unclear ?
                         in.indexmap().strides()[0], // stride of the in data
                         1,       // in : shift for multi fft.
                         out_fft, // out data
                         NULL,    // embed : unused. Doc unclear ?
                         out.indexmap().strides()[0], // stride of the out data
                         1, // out : shift for multi fft.
                         fftw_backward_forward, FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);
}

fourier_plan _fourier_base_plan(array_const_view<dcomplex, 2> in,
                                array_const_view<dcomplex, 2> out, int rank,
                                int *dims, int fftw_count,
                                int fftw_backward_forward) {

  auto p =
      fftw_plan_many_dft(rank,       // rank
                         dims,       // the dimension
                         fftw_count, // how many FFT : here 1
                         NULL,       // in data
                         NULL,       // embed : unused. Doc unclear ?
                         in.indexmap().strides()[0], // stride of the in data
                         1,    // in : shift for multi fft.
                         NULL, // out data
                         NULL, // embed : unused. Doc unclear ?
                         out.indexmap().strides()[0], // stride of the out data
                         1, // out : shift for multi fft.
                         fftw_backward_forward, FFTW_ESTIMATE);

  auto plan = std::make_unique<fourier_plan_base>((void *)p);

  return std::move(plan);
}

void _fourier_base(array_const_view<dcomplex, 2> in,
                   array_view<dcomplex, 2> out, fourier_plan &plan) {

  auto in_fft = reinterpret_cast<fftw_complex *>(in.data_start());
  auto out_fft = reinterpret_cast<fftw_complex *>(out.data_start());

  fftw_execute_dft((fftw_plan)plan.get()->plan_ptr, in_fft, out_fft);
}

void _fourier_base_destroy_plan(void *p) { fftw_destroy_plan((fftw_plan)p); }

fourier_plan_base::~fourier_plan_base() {
  _fourier_base_destroy_plan(plan_ptr);
}

} // namespace tprf::fourier
