/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, N. Wentzell, H. U.R. Strand
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

namespace tprf {
namespace fourier {

// template <template <typename, typename> class template_type> struct
// is_gf_view_or_cview : std::false_type {};
// template <> struct is_gf_view_or_cview<gf_view> : std::true_type {};
// template <> struct is_gf_view_or_cview<gf_const_view> : std::true_type {};

// 2D Fourier
template <typename M_out, typename M_in, typename T>
// First argument should be const &, but slicing has to be fixed first (creates
// const_view otherwise)
void _fourier_impl(
    gf_view<cartesian_product<M_out, M_out>, T> G_out,
    gf_const_view<cartesian_product<M_in, M_in>, T> const &G_in) {

  // Assume all in/out meshes to be equal and fetch the first
  auto mesh_in = std::get<0>(G_in.mesh());
  auto mesh_out = std::get<0>(G_out.mesh());

  TRIQS_ASSERT(G_out.target_shape() == G_in.target_shape());

  var_t _{};

  gf<cartesian_product<M_out, M_in>, T> G_1{{mesh_out, mesh_in},
                                            G_out.target_shape()};
  for (auto const &i : mesh_in)
    _fourier_impl(G_1[_][i], G_in[_][i]);
  for (auto const &o : mesh_out)
    _fourier_impl(G_out[o][_], make_const_view(G_1[o][_]));
}

// 3D Fourier
template <typename M_out, typename M_in, typename T>
// First argument should be const &, but slicing has to be fixed first (creates
// const_view otherwise)
void _fourier_impl(
    gf_view<cartesian_product<M_out, M_out, M_out>, T> G_out,
    gf_const_view<cartesian_product<M_in, M_in, M_in>, T> const &G_in) {

  // Assume all in/out meshes to be equal and fetch the first
  auto mesh_in = std::get<0>(G_in.mesh());
  auto mesh_out = std::get<0>(G_out.mesh());

  TRIQS_ASSERT(G_out.target_shape() == G_in.target_shape());

  var_t _{};

  gf<cartesian_product<M_out, M_in, M_in>, T> G_1{{mesh_out, mesh_in, mesh_in},
                                                  G_out.target_shape()};
  for (auto const &i1 : mesh_in)
    for (auto const &i2 : mesh_in)
      _fourier_impl(G_1[_][i1][i2], G_in[_][i1][i2]);

  gf<cartesian_product<M_out, M_out, M_in>, T> G_2{
      {mesh_out, mesh_out, mesh_in}, G_out.target_shape()};
  for (auto const &o : mesh_out)
    for (auto const &i : mesh_in)
      _fourier_impl(G_2[o][_][i], make_const_view(G_1[o][_][i]));

  for (auto const &o1 : mesh_out)
    for (auto const &o2 : mesh_out)
      _fourier_impl(G_out[o1][o2][_], make_const_view(G_2[o1][o2][_]));
}

// Block fourier
template <typename M_out, typename M_in, typename T>
void _fourier_impl(block_gf_view<M_out, T> const &G_out,
                   block_gf_const_view<M_in, T> const &G_in) {
  TRIQS_ASSERT(G_out.size() == G_in.size());
  for (int i : range(G_out.size()))
    _fourier_impl(G_out[i], G_in[i]);
}

// Block2 fourier
template <typename M_out, typename M_in, typename T>
void _fourier_impl(block2_gf_view<M_in, T> const &G_out,
                   block2_gf_const_view<M_out, T> const &G_in) {
  TRIQS_ASSERT(G_out.size1() == G_in.size1());
  TRIQS_ASSERT(G_out.size2() == G_in.size2());

  for (int i : range(G_out.size1()))
    for (int j : range(G_out.size2()))
      _fourier_impl(G_out(i, j), G_in(i, j));
}

template <template <typename, typename> class Gf, typename M_out, typename M_in,
          typename T>
typename Gf<M_out, T>::regular_type
_make_gf_impl(typename Gf<M_in, T>::const_view_type const &G_in, int n) {
  auto out_mesh = gf_mesh<M_out>{G_in.mesh().domain(), n};
  auto G_out = gf<M_out, T>{out_mesh, G_in.target_shape(), G_in.indices()};
  _fourier_impl(G_out(), G_in);
  return G_out;
}

template <template <typename, typename> class Gf, typename M_out, typename M_in,
          typename T>
typename Gf<cartesian_product<M_out, M_out>, T>::regular_type _make_gf_impl_2d(
    typename Gf<cartesian_product<M_in, M_in>, T>::const_view_type const &G_in,
    int n_1, int n_2) {
  auto domain = std::get<0>(G_in.mesh()).domain(); // Assume equal domains
  auto out_mesh =
      gf_mesh<cartesian_product<M_out, M_out>>{{domain, n_1}, {domain, n_2}};
  auto G_out = gf<cartesian_product<M_out, M_out>, T>{
      out_mesh, G_in.target_shape(), G_in.indices()};
  tprf::fourier::_fourier_impl(G_out(), G_in);
  return G_out;
}

template <template <typename, typename> class Gf, typename M_out, typename M_in,
          typename T>
typename Gf<cartesian_product<M_out, M_out, M_out>, T>::regular_type
_make_gf_impl_3d(typename Gf<cartesian_product<M_in, M_in, M_in>,
                             T>::const_view_type const &G_in,
                 int n_1, int n_2, int n_3) {
  auto domain = std::get<0>(G_in.mesh()).domain(); // Assume equal domains
  auto out_mesh = gf_mesh<cartesian_product<M_out, M_out, M_out>>{
      {domain, n_1}, {domain, n_2}, {domain, n_3}};
  auto G_out = gf<cartesian_product<M_out, M_out, M_out>, T>{
      out_mesh, G_in.target_shape(), G_in.indices()};
  tprf::fourier::_fourier_impl(G_out(), G_in);
  return G_out;
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<is_gf<Gf<imtime, T>>::value,
                 typename Gf<imfreq, T>::regular_type>
make_gf_from_fourier(Gf<imtime, T> const &G_tau, int n_iw = -1) {
  if (n_iw == -1)
    n_iw = (G_tau.mesh().size() - 1) / 2;
  return _make_gf_impl<Gf, imfreq, imfreq, T>(G_tau, n_iw);
}

// Factory to create imaginary time Green functions from inverse fourier
// transform of gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<is_gf<Gf<imfreq, T>>::value,
                 typename Gf<imtime, T>::regular_type>
make_gf_from_inverse_fourier(Gf<imfreq, T> const &G_iw, int n_tau = -1) {
  if (n_tau == -1)
    n_tau = 2 * (G_iw.mesh().last_index() + 1) + 1;
  return _make_gf_impl<Gf, imfreq, imtime, T>(G_iw, n_tau);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_gf<Gf<cartesian_product<imtime, imtime>, T>>::value,
    typename Gf<cartesian_product<imfreq, imfreq>, T>::regular_type>
make_gf_from_fourier(Gf<cartesian_product<imtime, imtime>, T> const &G_tau,
                     int n_iw1, int n_iw2) {
  return _make_gf_impl_2d<Gf, imfreq, imtime, T>(G_tau, n_iw1, n_iw2);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_gf<Gf<cartesian_product<imfreq, imfreq>, T>>::value,
    typename Gf<cartesian_product<imtime, imtime>, T>::regular_type>
make_gf_from_inverse_fourier(
    Gf<cartesian_product<imfreq, imfreq>, T> const &G_iw, int n_tau1,
    int n_tau2) {
  return _make_gf_impl_2d<Gf, imtime, imfreq, T>(G_iw, n_tau1, n_tau2);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_gf<Gf<cartesian_product<imtime, imtime, imtime>, T>>::value,
    typename Gf<cartesian_product<imfreq, imfreq, imfreq>, T>::regular_type>
make_gf_from_fourier(
    Gf<cartesian_product<imtime, imtime, imtime>, T> const &G_tau, int n_iw1,
    int n_iw2, int n_iw3) {
  return _make_gf_impl_3d<Gf, imfreq, imtime, T>(G_tau, n_iw1, n_iw2, n_iw3);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_gf<Gf<cartesian_product<imfreq, imfreq, imfreq>, T>>::value,
    typename Gf<cartesian_product<imtime, imtime, imtime>, T>::regular_type>
make_gf_from_inverse_fourier(
    Gf<cartesian_product<imfreq, imfreq, imfreq>, T> const &G_iw, int n_tau1,
    int n_tau2, int n_tau3) {
  return _make_gf_impl_3d<Gf, imtime, imfreq, T>(G_iw, n_tau1, n_tau2, n_tau3);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename M_out, typename M_in,
          typename T>
typename Gf<M_out, T>::regular_type
_make_block_gf_impl(typename Gf<M_in, T>::const_view_type const &G_in, int n) {
  // block_gf[_const][_view]
  if
    constexpr(Gf<M_in, T>::arity == 1) {
      std::vector<gf<M_out, T>> G_out_vec;
      // for (auto const &bl : G_in)
      // G_out_vec.emplace_back(make_gf_from_fourier(bl, n));
      for (auto const &bl : G_in)
        G_out_vec.emplace_back(
            _make_gf_impl<gf_const_view, M_out, M_in, T>(bl(), n));
      typename Gf<M_out, T>::regular_type G_out{G_in.block_names(), G_out_vec};
      return G_out;
    }
  // block2_gf[_const][_view]
  else if (Gf<M_in, T>::arity == 2) {
    std::vector<std::vector<gf<M_out, T>>> G_out_vecvec;
    for (int i : range(G_in.size1())) {
      std::vector<gf<M_out, T>> temp_vec;
      // for (int j : range(G_in.size2()))
      // temp_vec.emplace_back(make_gf_from_fourier(G_in(i,j), n));
      for (int j : range(G_in.size2()))
        temp_vec.emplace_back(
            _make_gf_impl<gf_const_view, M_out, M_in, T>(G_in(i, j), n));
      G_out_vecvec.emplace_back(std::move(temp_vec));
    }
    typename Gf<M_out, T>::regular_type G_out{G_in.block_names(), G_out_vecvec};
    return G_out;
  }
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of block[2]_gf_const_view
template <template <typename, typename> class Gf, typename M_out, typename M_in,
          typename T>
typename Gf<cartesian_product<M_out, M_out>, T>::regular_type
_make_block_gf_impl_2d(
    typename Gf<cartesian_product<M_in, M_in>, T>::const_view_type const &G_in,
    int n1, int n2) {

  using M_prod_in = cartesian_product<M_in, M_in>;
  using M_prod_out = cartesian_product<M_out, M_out>;

  // block_gf_const_view
  if
    constexpr(Gf<M_prod_in, T>::arity == 1) {
      std::vector<gf<M_prod_out, T>> G_out_vec;
      for (auto const &bl : G_in)
        G_out_vec.emplace_back(
            _make_gf_impl_2d<gf_const_view, M_out, M_in, T>(bl(), n1, n2));
      typename Gf<M_prod_out, T>::regular_type G_out{G_in.block_names(),
                                                     G_out_vec};
      return G_out;
    }
  // block2_gf_const_view
  else if (Gf<M_prod_in, T>::arity == 2) {
    std::vector<std::vector<gf<M_prod_out, T>>> G_out_vecvec;
    for (int i : range(G_in.size1())) {
      std::vector<gf<M_prod_out, T>> temp_vec;
      for (int j : range(G_in.size2()))
        temp_vec.emplace_back(_make_gf_impl_2d<gf_const_view, M_out, M_in, T>(
            G_in(i, j), n1, n2));
      G_out_vecvec.emplace_back(std::move(temp_vec));
    }
    typename Gf<M_prod_out, T>::regular_type G_out{G_in.block_names(),
                                                   G_out_vecvec};
    return G_out;
  }
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of block[2]_gf_const_view
template <template <typename, typename> class Gf, typename M_out, typename M_in,
          typename T>
typename Gf<cartesian_product<M_out, M_out, M_out>, T>::regular_type
_make_block_gf_impl_3d(typename Gf<cartesian_product<M_in, M_in, M_in>,
                                   T>::const_view_type const &G_in,
                       int n1, int n2, int n3) {

  using M_prod_in = cartesian_product<M_in, M_in, M_in>;
  using M_prod_out = cartesian_product<M_out, M_out, M_out>;

  // block_gf[_const][_view]
  if
    constexpr(Gf<M_prod_in, T>::arity == 1) {
      std::vector<gf<M_prod_out, T>> G_out_vec;
      for (auto const &bl : G_in)
        G_out_vec.emplace_back(
            _make_gf_impl_3d<gf_const_view, M_out, M_in, T>(bl(), n1, n2, n3));
      typename Gf<M_prod_out, T>::regular_type G_out{G_in.block_names(),
                                                     G_out_vec};
      return G_out;
    }
  // block2_gf[_const][_view]
  else if (Gf<M_prod_in, T>::arity == 2) {
    std::vector<std::vector<gf<M_prod_out, T>>> G_out_vecvec;
    for (int i : range(G_in.size1())) {
      std::vector<gf<M_prod_out, T>> temp_vec;
      for (int j : range(G_in.size2()))
        temp_vec.emplace_back(_make_gf_impl_3d<gf_const_view, M_out, M_in, T>(
            G_in(i, j), n1, n2, n3));
      G_out_vecvec.emplace_back(std::move(temp_vec));
    }
    typename Gf<M_prod_out, T>::regular_type G_out{G_in.block_names(),
                                                   G_out_vecvec};
    return G_out;
  }
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<is_block_gf_or_view<Gf<imtime, T>>::value,
                 typename Gf<imfreq, T>::regular_type>
make_gf_from_fourier(Gf<imtime, T> const &G_tau, int n_iw = -1) {
  return _make_block_gf_impl<Gf, imfreq, imtime, T>(G_tau, n_iw);
}

// Factory to create imaginary time Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<is_block_gf_or_view<Gf<imfreq, T>>::value,
                 typename Gf<imtime, T>::regular_type>
make_gf_from_inverse_fourier(Gf<imfreq, T> const &G_iw, int n_tau = -1) {
  return _make_block_gf_impl<Gf, imtime, imfreq, T>(G_iw, n_tau);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_block_gf_or_view<Gf<cartesian_product<imtime, imtime>, T>>::value,
    typename Gf<cartesian_product<imfreq, imfreq>, T>::regular_type>
make_gf_from_fourier(Gf<cartesian_product<imtime, imtime>, T> const &G_tau,
                     int n_iw1, int n_iw2) {
  return _make_block_gf_impl_2d<Gf, imfreq, imtime, T>(G_tau, n_iw1, n_iw2);
}

// Factory to create imaginary time Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_block_gf_or_view<Gf<cartesian_product<imfreq, imfreq>, T>>::value,
    typename Gf<cartesian_product<imtime, imtime>, T>::regular_type>
make_gf_from_inverse_fourier(
    Gf<cartesian_product<imfreq, imfreq>, T> const &G_iw, int n_tau1,
    int n_tau2) {
  return _make_block_gf_impl_2d<Gf, imtime, imfreq, T>(G_iw, n_tau1, n_tau2);
}

// Factory to create Matsubara frequency Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_block_gf_or_view<
        Gf<cartesian_product<imtime, imtime, imtime>, T>>::value,
    typename Gf<cartesian_product<imfreq, imfreq, imfreq>, T>::regular_type>
make_gf_from_fourier(
    Gf<cartesian_product<imtime, imtime, imtime>, T> const &G_tau, int n_iw1,
    int n_iw2, int n_iw3) {
  return _make_block_gf_impl_3d<Gf, imfreq, imtime, T>(G_tau, n_iw1, n_iw2,
                                                       n_iw3);
}

// Factory to create imaginary time Green functions from inverse fourier
// transform of block[2]_gf[_const][_view]
template <template <typename, typename> class Gf, typename T>
std::enable_if_t<
    is_block_gf_or_view<Gf<cartesian_product<imfreq, imfreq>, T>>::value,
    typename Gf<cartesian_product<imtime, imtime, imtime>, T>::regular_type>
make_gf_from_inverse_fourier(
    Gf<cartesian_product<imfreq, imfreq, imfreq>, T> const &G_iw, int n_tau1,
    int n_tau2, int n_tau3) {
  return _make_block_gf_impl_3d<Gf, imtime, imfreq, T>(G_iw, n_tau1, n_tau2,
                                                       n_tau3);
}

} // namespace fourier
} // namespace tprf
