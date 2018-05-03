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
#pragma once

#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>

using namespace triqs::gfs;
using namespace triqs::arrays;

#include "types.hpp"

namespace tprf {

// ----------------------------------------------------

template <Channel_t C> class channel_grouping {
public:
  Channel_t channel = C;

  inline memory_layout_t<6> memory_layout() const;
  inline matrix_view<scalar_t> matrix_view(array_view<scalar_t, 6> arr) const;
};

// ----------------------------------------------------
// Index grouping definitions
// ----------------------------------------------------

// The standard layout is:
//
// {nu_1, nu_2, a, b, c, d} <=> {0, 1, 2, 3, 4, 5}
//
// where nu_1 and nu_2 (nu_1=0, nu_2=1) are fermionic Matsubara frequencies
// and a, b, c, d (a=2, b=3, c=4, d=5) are target-space indices

// ----------------------------------------------------
// Channel_t::PH

// in the particle-hole channel (Channel_t::PH) the indices are grouped as
// {nu_1, a, b}, {nu_2, c, d} <=> {0, 2, 3}, {1, 4, 5}

template <>
inline memory_layout_t<6>
channel_grouping<Channel_t::PH>::memory_layout() const {
  return make_memory_layout(0, 2, 3, 1, 4, 5);
}

template <>
inline matrix_view<scalar_t> channel_grouping<Channel_t::PH>::matrix_view(
    array_view<scalar_t, 6> arr) const {
  return make_matrix_view(group_indices_view(arr, {0, 2, 3}, {1, 4, 5}));
}

// ----------------------------------------------------
// Channel_t::PH_bar

// in the particle-hole-bar channel (Channel_t::PH_bar) the indices are grouped
// as
// {nu_1, a, d}, {nu_2, c, b} <=> {0, 2, 5}, {1, 4, 3}

template <>
inline memory_layout_t<6>
channel_grouping<Channel_t::PH_bar>::memory_layout() const {
  return make_memory_layout(0, 2, 5, 1, 4, 3);
}

template <>
inline matrix_view<scalar_t> channel_grouping<Channel_t::PH_bar>::matrix_view(
    array_view<scalar_t, 6> arr) const {
  return make_matrix_view(group_indices_view(arr, {0, 2, 5}, {1, 4, 3}));
}

// ----------------------------------------------------
// Channel_t::PP

// in the particle-particle channel (Channel_t::PP) the indices are grouped as
// {nu_1, a, c}, {nu_2, b, d} <=> {0, 2, 4}, {1, 3, 5}

template <>
inline memory_layout_t<6>
channel_grouping<Channel_t::PP>::memory_layout() const {
  return make_memory_layout(0, 2, 4, 1, 3, 5);
}

template <>
inline matrix_view<scalar_t> channel_grouping<Channel_t::PP>::matrix_view(
    array_view<scalar_t, 6> arr) const {
  return make_matrix_view(group_indices_view(arr, {0, 2, 4}, {1, 3, 5}));
}

} // namespace tprf
