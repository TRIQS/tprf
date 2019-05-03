/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018, The Simons Foundation
 * Author: H. U.R. Strand
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
#pragma once

#include <itertools/itertools.hpp>
#include <mpi/mpi.hpp>

#include "types.hpp"

namespace triqs_tprf {

template<class T>
auto mpi_view(const array<T, 1> &arr, mpi::communicator const & c) {

  auto slice = itertools::chunk_range(0, arr.shape()[0], c.size(), c.rank());

  /*
  std::cout << "mpi_view<array> " << "rank = " << c.rank() << " size = " << arr.shape()[0]
  	    << " s,e = " << slice.first << ", " << slice.second << "\n";
  */

  return arr(range(slice.first, slice.second + 1));
}

template<class T>
auto mpi_view(const array<T, 1> &arr) {
  mpi::communicator c;
  return mpi_view(arr, c);
}
  
template<class T>
auto mpi_view(const gf_mesh<T> &mesh, mpi::communicator const & c) {

  auto slice = itertools::chunk_range(0, mesh.size(), c.size(), c.rank());
  int size = slice.second - slice.first;

  /*
  std::cout << "mpi_view<mesh> " << "rank = " << c.rank() << " size = " << size
  	    << " s,e = " << slice.first << ", " << slice.second << "\n";
  c.barrier();
  */
  
  array<typename gf_mesh<T>::mesh_point_t, 1> arr(size);

  auto iter = mesh.begin();
  iter += slice.first;

  for ( auto idx : range(0, size) ) {
    auto w = *iter;
    arr(idx) = w;
    iter++;
  }

  /*
  std::cout << "mpi_view<mesh> " << "rank = " << c.rank() << " size = " << size
  	    << " s,e = " << slice.first << ", " << slice.second << " arr size " << arr.size() << "\n";
  c.barrier();
  */

  return arr;
}

template<class T>
auto mpi_view(const gf_mesh<T> &mesh) {
  mpi::communicator c;
  return mpi_view(mesh, c);
}
  
} // namespace triqs_tprf
