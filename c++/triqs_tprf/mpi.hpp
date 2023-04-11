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

#include <span>

#include "types.hpp"

namespace triqs_tprf {

  template <class T> auto mpi_view(const T &mesh, mpi::communicator c = {}) {

    auto slice = itertools::chunk_range(0, mesh.size(), c.size(), c.rank());
    int size   = slice.second - slice.first;

    std::vector<typename T::mesh_point_t> arr{};
    arr.reserve(size);

    auto iter = std::next(mesh.begin(), slice.first);
    for (auto idx : range(0, size)) {
      arr.emplace_back(*iter);
      iter++;
    }

    return arr;
  }

} // namespace triqs_tprf
