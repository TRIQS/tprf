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

namespace triqs_tprf::fourier {

  /*
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
  */

fourier_plan _fourier_base_plan(array_const_view<dcomplex, 2> in,
                                array_const_view<dcomplex, 2> out, int rank,
                                int *dims, int fftw_count,
                                int fftw_backward_forward) {

  auto in_fft  = reinterpret_cast<fftw_complex *>(in.data_start());
  auto out_fft = reinterpret_cast<fftw_complex *>(out.data_start());

  /*
  std::cout << "--> triqs_tprf::fourier _fourier_base_plan\n";
  std::cout << "rank = " << rank << "\n";
  std::cout << "dims[0] = " << dims[0] << "\n";
  std::cout << "dims[1] = " << dims[1] << "\n";
  std::cout << "dims[2] = " << dims[2] << "\n";
  std::cout << "fftw_count = " << fftw_count << "\n";
  std::cout << "in_strides = " << in.indexmap().strides() << "\n";
  std::cout << "out_strides = " << out.indexmap().strides() << "\n";
  */
  
  auto p =
      fftw_plan_many_dft(rank,       // rank
                         dims,       // the dimension
                         fftw_count, // how many FFT : here 1
                         //NULL,       // in data
			 in_fft,
                         NULL,       // embed : unused. Doc unclear ?
                         in.indexmap().strides()[0], // stride of the in data
                         1,    // in : shift for multi fft.
                         //NULL, // out data
			 out_fft,
                         NULL, // embed : unused. Doc unclear ?
                         out.indexmap().strides()[0], // stride of the out data
                         1, // out : shift for multi fft.
                         fftw_backward_forward, FFTW_ESTIMATE);

  //fftw_print_plan(p); std::cout << "\n";
  auto plan = std::make_unique<fourier_plan_base>((void *)p);
  //fftw_print_plan((fftw_plan)plan.get()->plan_ptr); std::cout << "\n";

  return std::move(plan);
}

void _fourier_base(array_const_view<dcomplex, 2> in,
                   array_view<dcomplex, 2> out, fourier_plan &plan) {

  auto in_fft = reinterpret_cast<fftw_complex *>(in.data_start());
  auto out_fft = reinterpret_cast<fftw_complex *>(out.data_start());

  //std::cout << "in = " << in << "\n";
  //std::cout << "out = " << out << "\n";

  //fftw_print_plan((fftw_plan)plan.get()->plan_ptr); std::cout << "\n";
  
  fftw_execute_dft((fftw_plan)plan.get()->plan_ptr, in_fft, out_fft);

  //std::cout << "out = " << out << "\n";
  
}

void _fourier_base_destroy_plan(void *p) { fftw_destroy_plan((fftw_plan)p); }

fourier_plan_base::~fourier_plan_base() {
  _fourier_base_destroy_plan(plan_ptr);
}

} // namespace triqs_tprf::fourier
