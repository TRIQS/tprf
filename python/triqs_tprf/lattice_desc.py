# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/lattice.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m lattice -o lattice -C pytriqs --moduledoc="Lattice functionality" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "lattice", doc = r"Lattice functionality", app_name = "triqs_tprf")

# Imports
module.add_imports(*['pytriqs.gf', 'pytriqs.lattice'])

# Add here all includes
module.add_include("triqs_tprf/lattice.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/tuple.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", doc = r"""Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g0_wk (double mu, triqs_tprf::e_k_cvt e_k, gf_mesh<triqs::gfs::imfreq> mesh)", doc = r"""Construct a non-interacting Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`

  Computes

  .. math::
     G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
         (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
	 \right]^{-1}_{a\bar{b}},

  using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, chemical potential :math:`\mu`,
  and a Matsubara frequency Green's function mesh.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

mesh
     imaginary frequency mesh

Returns
-------
out
     Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_w_cvt sigma_w)", doc = r"""Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a momentum independent Matsubara frequency
 self energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

sigma_w
     imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`

Returns
-------
out
     Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_wk_cvt sigma_wk)", doc = r"""Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a momentum independent Matsubara frequency
 self energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

sigma_wk
     imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`

Returns
-------
out
     Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_w_t triqs_tprf::lattice_dyson_g_w (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_w_cvt sigma_w)", doc = r"""Construct an interacting Matsubara frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function :math:`G_{a\bar{b}}(i\omega_n)`

 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n) = \frac{1}{N_k} \sum_\mathbf{k} \left[
        (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
	\right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a momentum independent Matsubara frequency
 self energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

sigma_w
     imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`

Returns
-------
out
     Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wr_t triqs_tprf::fourier_wk_to_wr (triqs_tprf::g_wk_cvt g_wk)", doc = r"""Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{k}) \right\}`

Parameters
----------
g_wk
     k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

Returns
-------
out
     real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::fourier_wr_to_wk (triqs_tprf::g_wr_cvt g_wr)", doc = r"""Fast fourier transform of imaginary frequency Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

Parameters
----------
g_wr
     real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`

Returns
-------
out
     k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_tr_t triqs_tprf::fourier_wr_to_tr (triqs_tprf::g_wr_cvt g_wr, int nt = -1)", doc = r"""Fast fourier transform of real-space Green's function from Matsubara frequency to imaginary time

    Computes: :math:`G_{a\bar{b}}(\tau, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

Parameters
----------
g_wr
     real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`

Returns
-------
out
     real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`""")

module.add_function ("triqs_tprf::g_wr_t triqs_tprf::fourier_tr_to_wr (triqs_tprf::g_tr_cvt g_tr, int nw = -1)", doc = r"""Fast fourier transform of real-space Green's function from imaginary time to Matsubara frequency

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(\tau, \mathbf{r}) \right\}`

Parameters
----------
g_tr
     real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`

Returns
-------
out
     real-space Matsubara frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::lindhard_chi00_wk (triqs_tprf::e_k_cvt e_k, int nw, double beta, double mu)", doc = r"""Generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`.

   Analytic calculation of the generalized (non-interacting) Lindhard susceptibility
   in the particle-hole channel. The analytic expression is obtained using residue calculus
   to explicitly evaluate the matsubara sum of the fourier transformed imaginary time
   bubble product of two non-interacting single-particle Green's functions.

   .. math::
      G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) =
      \left[ i\omega_n \cdot \mathbf{1} - \epsilon(\mathbf{k}) \right]^{-1} .

   The analytic evaluation of the bubble diagram gives

   .. math::
        \chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q}) \equiv
        \mathcal{F} \left\{
	  - G^{(0)}_{d\bar{a}}(\tau, \mathbf{r}) G^{(0)}_{b\bar{c}}(-\tau, -\mathbf{r})
	\right\}
        =
	- \frac{1}{N_k} \sum_{\nu} \sum_{\mathbf{k}}
          G^{(0)}_{d\bar{a}}(\nu, \mathbf{k})
	  G^{(0)}_{b\bar{c}}(\nu + \omega, \mathbf{k} + \mathbf{q})
	\\ =
	- \frac{1}{N_k} \sum_{\nu} \sum_{\mathbf{k}}
	  \left( \sum_{i}
          U^\dagger_{di}(\mathbf{k}) \frac{1}{i\nu - \epsilon_{\mathbf{k}, i}} U_{i\bar{a}}(\mathbf{k})
	  \right)
	  \left( \sum_j
	  U^\dagger_{bj}(\mathbf{k} + \mathbf{q})
	  \frac{1}{i\nu + i\omega - \epsilon_{\mathbf{k} + \mathbf{q}, j}}
	  U_{j\bar{c}}(\mathbf{k} + \mathbf{q})
	  \right)
	\\ =
	\frac{1}{N_k} \sum_{\mathbf{k}} \sum_{ij}
	  \left(
	    [1 - \delta_{0, \omega_n} \delta_{\epsilon_{\mathbf{k},i},\epsilon_{\mathbf{k}+\mathbf{q}, j}})]
	    \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
	         {i\omega_n + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
	    +
  	    \delta_{0, \omega_n} \delta_{\epsilon_{\mathbf{k},i},\epsilon_{\mathbf{k}+\mathbf{q}, j}}
	    \frac{\beta}{4 \cosh^2 (\beta \epsilon_{\mathbf{k}, i} / 2) }
	  \right)
	  \\ \times
	  U_{i\bar{a}}(\mathbf{k}) U^\dagger_{di}(\mathbf{k})
	  U_{j\bar{c}}(\mathbf{k} + \mathbf{q}) U^\dagger_{bj}(\mathbf{k} + \mathbf{q})

   where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   .. math::
      \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
      = \delta_{ij} \epsilon_{\mathbf{k}, i}

   .. note::
      The analytic formula is sub-optimal in terms of performance for higher temperatures. The evaluation
      scales as :math:`\mathcal{O}(N_k^2)` which is worse than computing the bubble explicitly in imaginary
      time, with scaling :math:`\mathcal{O}(N_k N_\tau \log(N_k N_\tau)` for :math:`N_k \gg N_\tau`.

   .. note::
      Care must be taken when evaluating the fermionic Matsubara frequency sum of the
      product of two simple poles. By extending the sum to an integral over the complex
      plane the standard expression for the Lindhard response is obtained when the
      poles are non-degenerate. The degenerate case produces an additional frequency independent
      contribution (the last term on the last row).""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::solve_rpa_PH (triqs_tprf::chi_wk_vt chi0, array_view<std::complex<double>,4> U)", doc = r"""Random Phase Approximation (RPA) in the particle-hole channel

     Computes the equation

     .. math::
         \chi(\bar{a}b\bar{c}d) = \big(
         \mathbb{1}
         - \chi^{(0)}(\bar{a}b\bar{B}A) U(A\bar{B}D\bar{C})
         \big)^{-1} \chi^{(0)}(\bar{C}D\bar{c}d)\,.

Parameters
----------
chi0
     bare particle-hole bubble :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

U
     RPA static vertex as obtained from triqs_tprf.rpa_tensor.get_rpa_tensor :math:`U_{a\bar{b}c\bar{d}}`

Returns
-------
out
     RPA suceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::dynamical_screened_interaction_W_wk (triqs_tprf::chi_wk_cvt PI_wk, triqs_tprf::chi_k_cvt V_k)", doc = r"""Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator
    for static momentum-dependent interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

    Instead of returning :math:`W^{(full)}` we return the dynamical/retarded part
    :math:`W^{(r)}` (with zero high-frequency offset)

    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) =
            W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) - V_{abcd}(\mathbf{k})

Parameters
----------
PI_wk
     polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

V_k
     static interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::dynamical_screened_interaction_W_wk_from_generalized_susceptibility (triqs_tprf::chi_wk_cvt chi_wk, triqs_tprf::chi_k_cvt V_k)", doc = r"""Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator
    for static momentum-dependent interactions :math:`V(\mathbf{k})` and
    known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
          V_{abcd}(\mathbf{k}) +
	  \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
          \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
          V_{hgcd}(\mathbf{k})

    Instead of returning :math:`W^{(full)}` we return the dynamical/retarded part
    :math:`W^{(r)}` (with zero high-frequency offset)

    .. math::
        W_{abcd}(i\omega_n, \mathbf{k}) =
            W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) - V_{abcd}(\mathbf{k})

Parameters
----------
chi_wk
     polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

V_k
     static interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::gw_sigma_wk_serial_fft (triqs_tprf::chi_wk_cvt Wr_wk, triqs_tprf::g_wk_cvt g_wk)", doc = r"""GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator

    Fourier transforms the screened interaction and the single-particle
    Green's function to imagiary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W^{(r)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{(r)}_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W^{(r)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma_{ab}(\tau, \mathbf{r}) \right\}

Parameters
----------
V_k
     static bare interaction :math:`V_{abcd}(\mathbf{k})`

Wr_wk
     retarded screened interaction :math:`W^{(r)}_{abcd}(i\omega_n, \mathbf{k})`

g_wk
     single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
out
     GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_tr_t triqs_tprf::gw_sigma_tr (triqs_tprf::chi_tr_cvt Wr_tr, triqs_tprf::g_tr_cvt g_tr)", doc = r"""GW self energy :math:`\Sigma(\tau, \mathbf{r})` calculator

    Computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          \sum_{cd} W^{(r)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

Parameters
----------
Wr_tr
     retarded screened interaction :math:`W^{(r)}_{abcd}(\tau, \mathbf{r})`

g_tr
     single particle Green's function :math:`G_{ab}(\tau, \mathbf{r})`

Returns
-------
out
     GW self-energy :math:`\Sigma_{ab}(\tau, \mathbf{r})`""")

module.add_function ("triqs_tprf::gk_iw_t triqs_tprf::eliashberg_product (triqs_tprf::chi_wk_vt Gamma_pp, triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)", doc = r"""Linearized Eliashberg product

     Computes the product

     .. math::
         \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{k}'} \sum_{i\nu'}
	 \Gamma_{A\bar{a}B\bar{b}}(\mathbf{k}-\mathbf{k}', i\nu - i\nu')
	 \\ \times
	 G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')

Parameters
----------
chi_pp
     particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\nu_n)`

g_kw
     single particle Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\nu_n)`

delta_kw
     pairing self-energy :math:`\Delta_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n)`

Returns
-------
out
     Gives the result of the product :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG \Delta`""")

module.add_function ("triqs_tprf::gk_iw_t triqs_tprf::eliashberg_product_fft (triqs_tprf::chi_tr_vt Gamma_pp_dyn_tr, triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)", doc = r"""Linearized Eliashberg product via FFT

     Computes the product

     .. math::
         \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =  -\frac{1}{N_k \beta}\sum_{\mathbf{k}'} \sum_{i\nu'}
	 \Gamma_{A\bar{a}B\bar{b}}(\mathbf{k}-\mathbf{k}', i\nu - i\nu')
	 \\ \times
	 G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')\,,

     by taking advantage of the convolution theorem.

     We therefore first calculate

     .. math::
        \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{r}, \tau) =
	 -\Gamma_{A\bar{a}B\bar{b}}(\mathbf{r}, \tau) F_{AB}(\mathbf{r}, \tau) \,,

     where

     .. math::
        F_{AB}(\mathbf{r}, \tau)  =
        \mathcal{F}\big(G_{A\bar{c}}(\mathbf{k}', i\nu')
	 \Delta_{\bar{c}\bar{d}}(\mathbf{k}', i\nu')
	 G_{B\bar{d}}(-\mathbf{k}', -i\nu')\big)\,.

     Then we Fourier transform

     .. math::
          \Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{k},i\nu) =
          \mathcal{F}\big(\Delta^{(out)}_{\bar{a}\bar{b}}(\mathbf{r}, \tau)\big)\,,

    to get the same result, but with far less computational effort.

Parameters
----------
chi_rt
     dynamic part of the particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{r}, \tau)`

chi_r
     constant part of the particle-particle vertex :math:`\Gamma^{(pp)}_{a\bar{b}c\bar{d}}(\mathbf{r})`

g_kw
     single particle Green's function :math:`G_{a\bar{b}}(\mathbf{k}, i\nu_n)`

delta_kw
     pairing self-energy :math:`\Delta_{\bar{a}\bar{b}}(\mathbf{k}, i\nu_n)`

Returns
-------
out
     Gives the result of the product :math:`\Delta^{(out)} \sim \Gamma^{(pp)}GG \Delta`""")

module.add_function ("triqs_tprf::gk_iw_t triqs_tprf::eliashberg_g_delta_g_product (triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)", doc = r"""""")

module.add_function ("std::tuple<chi_wk_vt,chi_k_vt> triqs_tprf::split_into_dynamic_wk_and_constant_k (triqs_tprf::chi_wk_vt Gamma_pp)", doc = r"""""")

module.add_function ("std::tuple<chi_tr_vt,chi_r_vt> triqs_tprf::dynamic_and_constant_to_tr (triqs_tprf::chi_wk_vt Gamma_pp_dyn_wk, triqs_tprf::chi_k_vt Gamma_pp_const_k)", doc = r"""""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::gamma_PP_singlet (triqs_tprf::chi_wk_vt chi_c, triqs_tprf::chi_wk_vt chi_s, array_view<std::complex<double>,4> U_c, array_view<std::complex<double>,4> U_s)", doc = r"""Gamma particle-particle singlet

     Computes the particle-particle vertex for singlet pairing in the RPA limit

    .. math::
        \Gamma^{(\mathrm{singlet})}(a\bar{b}c\bar{d}) =
        \frac{3}{2} U^{(\mathrm{s})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{s})}(\bar{B}A\bar{C}D)
        U^{(\mathrm{s})}(D\bar{C}c\bar{d}) \\
        -\frac{1}{2} U^{(\mathrm{c})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{c})}(\bar{B}A\bar{C}D)
        U^{(\mathrm{c})}(D\bar{C}c\bar{d}) \\
       + \frac{1}{2}\big(U^{(\mathrm{s})}(a\bar{b}c\bar{d})+
        U^{(\mathrm{c})}(a\bar{b}c\bar{d})\big)

Parameters
----------
chi_c
     charge susceptibility  :math:`\chi^{(\mathrm{c})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

chi_s
     spin susceptibility  :math:`\chi^{(\mathrm{s})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

U_c
     charge interaction  :math:`U^{(\mathrm{c})}_{a\bar{b}c\bar{d}}`

U_s
     spin interaction  :math:`U^{(\mathrm{s})}_{a\bar{b}c\bar{d}}`

Returns
-------
out
     :math:`\Gamma^{(\mathrm{singlet})}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\omega_n)`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::gamma_PP_triplet (triqs_tprf::chi_wk_vt chi_c, triqs_tprf::chi_wk_vt chi_s, array_view<std::complex<double>,4> U_c, array_view<std::complex<double>,4> U_s)", doc = r"""Gamma particle-particle triplet

     Computes the particle-particle vertex for triplet pairing in the RPA limit

    .. math::
        \Gamma^{(\mathrm{triplet})}(a\bar{b}c\bar{d}) =
        -\frac{1}{2} U^{(\mathrm{s})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{s})}(\bar{B}A\bar{C}D)
        U^{(\mathrm{s})}(D\bar{C}c\bar{d}) \\
        -\frac{1}{2} U^{(\mathrm{c})}(a\bar{b}A\bar{B}) \chi^{(\mathrm{c})}(\bar{B}A\bar{C}D)
        U^{(\mathrm{c})}(D\bar{C}c\bar{d}) \\
       + \frac{1}{2}\big(U^{(\mathrm{s})}(a\bar{b}c\bar{d})+
        U^{(\mathrm{c})}(a\bar{b}c\bar{d})\big)

Parameters
----------
chi_c
     charge susceptibility  :math:`\chi^{(\mathrm{c})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

chi_s
     spin susceptibility  :math:`\chi^{(\mathrm{s})}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

U_c
     charge interaction  :math:`U^{(\mathrm{c})}_{a\bar{b}c\bar{d}}`

U_s
     spin interaction  :math:`U^{(\mathrm{s})}_{a\bar{b}c\bar{d}}`

Returns
-------
out
     :math:`\Gamma^{(\mathrm{triplet})}_{a\bar{b}c\bar{d}}(\mathbf{k}, i\omega_n)`""")

module.add_function ("array<std::complex<double>,6> triqs_tprf::cluster_mesh_fourier_interpolation (array<double,2> k_vecs, triqs_tprf::chi_wr_cvt chi)", doc = r"""""")

module.add_function ("triqs_tprf::chi_tr_t triqs_tprf::chi0_tr_from_grt_PH (triqs_tprf::g_tr_cvt g_tr)", doc = r"""Generalized susceptibility imaginary time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
     - G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})

Parameters
----------
g_tr
     Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi0_w0r_from_grt_PH (triqs_tprf::g_tr_cvt g_tr)", doc = r"""Generalized susceptibility zero imaginary frequency bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r}) =
     - \int_0^\beta d\tau \,
     G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})

Parameters
----------
g_tr
     Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r})` in real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi_w0r_from_chi_tr (triqs_tprf::chi_tr_cvt chi_tr)", doc = r"""Static susceptibility calculation :math:`\chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`

  Explicit calculation of the static, zero frequency response, by 2nd order trapetzoidal
  integration in imaginary time, i.e.

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) =
         \int_0^\beta d\tau \, \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})

Parameters
----------
chi_tr
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`
     in imaginary time and real space.

Returns
-------
out
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`
     at zero Matsubara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi_wr_from_chi_tr (triqs_tprf::chi_tr_cvt chi_tr, int nw)", doc = r"""Parallell Fourier transform from  :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \mathcal{F}_{\tau \rightarrow \omega} \left\{
         \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})
         \right\}

Parameters
----------
chi_tr
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`
     in imaginary time and real space.

Returns
-------
out
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
     in Matsubara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_tr_t triqs_tprf::chi_tr_from_chi_wr (triqs_tprf::chi_wr_cvt chi_wr, int ntau = -1)", doc = r"""Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`

  Computes

  .. math::
         \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
         \mathcal{F}_{\omega \rightarrow \tau} \left\{
	 \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \right\}

Parameters
----------
chi_tr
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`
     in imaginary time and real space.

Returns
-------
out
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
     in Matsubara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::chi_wk_from_chi_wr (triqs_tprf::chi_wr_cvt chi_wr)", doc = r"""Parallell Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
         \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{k}} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})
         \right\}

Parameters
----------
chi_wr
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
     in Matsubara frequency and real space.

Returns
-------
out
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`
     in Matsubara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi_wr_from_chi_wk (triqs_tprf::chi_wk_cvt chi_wk)", doc = r"""Parallell Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
         \mathcal{F}_{\mathbf{k} \rightarrow \mathbf{r}} \left\{
         \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})
         \right\}

Parameters
----------
chi_wr
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`
     in imaginary time and momentum space.

Returns
-------
out
     Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`
     in Matsubara frequency and real space.""")

module.add_function ("chi_t_t::zero_t triqs_tprf::chi_trapz_tau (triqs_tprf::chi_t_cvt chi_t)", doc = r"""""")

module.add_function ("triqs_tprf::chi_wnr_t triqs_tprf::chi0r_from_gr_PH (int nw, int nn, triqs_tprf::g_wr_cvt g_nr)", doc = r"""Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})`.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})

Parameters
----------
nw
     Number of bosonic Matsubara freqiencies.

nn
     Number of fermionic Matsubara freqiencies.

g_tr
     Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wnr_t triqs_tprf::chi0r_from_gr_PH_nompi (int nw, int nn, triqs_tprf::g_wr_cvt g_nr)", doc = r"""Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` without MPI parallellization.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})

Parameters
----------
nw
     Number of bosonic Matsubara freqiencies.

nn
     Number of fermionic Matsubara freqiencies.

g_tr
     Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wnk_t triqs_tprf::chi0q_from_g_wk_PH (int nw, int nn, triqs_tprf::g_wk_cvt g_wk)", doc = r"""Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` with convolution in k-space.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
     - \frac{\beta}{N_k} \sum_\mathbf{k}
     G_{d\bar{a}}(\nu, \mathbf{k}) \cdot G_{b\bar{c}}(\nu + \omega, \mathbf{k} - \mathbf{q})

Parameters
----------
nw
     Number of bosonic Matsubara freqiencies.

nn
     Number of fermionic Matsubara freqiencies.

g_tr
     Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wnr_t triqs_tprf::chi0r_from_chi0q (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = r"""Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum-space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real-space.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
     \mathcal{F}_{\mathbf{q} \rightarrow \mathbf{r}} \left\{
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})
     \right\}

Parameters
----------
chi_wnk
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.""")

module.add_function ("triqs_tprf::chi_wnk_t triqs_tprf::chi0q_from_chi0r (triqs_tprf::chi_wnr_cvt chi_wnr)", doc = r"""Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum space.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
     \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{q}} \left\{
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})
     \right\}

Parameters
----------
chi_wnr
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::chi0q_sum_nu (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = r"""Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max} \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

Parameters
----------
chi_wnk
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
out
     Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::chi0q_sum_nu_tail_corr_PH (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = r"""Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})` using higher order tail corrections when summing to infinity.

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{\beta^2} \sum_{\nu=-\infty}^\infty \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

Parameters
----------
chi_wnk
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
out
     Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_w_t triqs_tprf::chi0q_sum_nu_q (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = r"""Sum over fermionic frequency and momentum in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \frac{1}{N_k} \sum_\matbf{k} \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max}
     \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

Parameters
----------
chi_wnk
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
out
     Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega)` in one bosonic Matsubara frequency.""")

module.add_function ("triqs_tprf::chi_kwnn_t triqs_tprf::chiq_from_chi0q_and_gamma_PH (triqs_tprf::chi_wnk_cvt chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn)", doc = r"""Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k})`.

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k}) =
     \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}

Parameters
----------
chi0_wnk
     Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

gamma_ph_wnn
     Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.

Returns
-------
out
     Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k})`.""")

module.add_function ("triqs_tprf::chi_kw_t triqs_tprf::chiq_sum_nu_from_chi0q_and_gamma_PH (triqs_tprf::chi_wnk_cvt chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn)", doc = r"""Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

  Computes

  .. math::
     \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
     \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}

Parameters
----------
chi0_wnk
     Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

gamma_ph_wnn
     Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.

Returns
-------
out
     Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.""")

module.add_function ("gf<cartesian_product<triqs::lattice::brillouin_zone,triqs::gfs::imfreq>,tensor_valued<4>> triqs_tprf::chiq_sum_nu_from_g_wk_and_gamma_PH (triqs_tprf::gk_iw_t g_wk, triqs_tprf::g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1)", doc = r"""""")

module.add_function ("gf<cartesian_product<triqs::lattice::brillouin_zone,triqs::gfs::imfreq>,tensor_valued<4>> triqs_tprf::chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH (double mu, triqs_tprf::ek_vt e_k, triqs_tprf::g_iw_vt sigma_w, triqs_tprf::g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1)", doc = r"""""")

module.add_function ("gf<cartesian_product<triqs::lattice::brillouin_zone,triqs::gfs::imfreq>,tensor_valued<4>> triqs_tprf::chiq_sum_nu (triqs_tprf::chiq_t chiq)", doc = r"""""")

module.add_function ("gf<triqs::gfs::imfreq,tensor_valued<4>> triqs_tprf::chiq_sum_nu_q (triqs_tprf::chiq_t chiq)", doc = r"""""")



module.generate_code()