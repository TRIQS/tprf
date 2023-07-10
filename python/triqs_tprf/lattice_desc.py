# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/lattice.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m lattice -o lattice -C triqs -C nda_py --includes=../../c++ --moduledoc="Lattice functionality" --cxxflags="-std=c++20"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "lattice", doc = r"Lattice functionality", app_name = "triqs_tprf")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes'])

# Add here all includes
module.add_include("triqs_tprf/lattice.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/tuple.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/mesh.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", doc = r"""Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g0_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs::mesh::imfreq mesh)", doc = r"""Construct a non-interacting Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`

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

module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::lattice_dyson_g0_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs::mesh::dlr_imfreq mesh)")

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

module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_Dw_cvt sigma_w)")
                     
module.add_function ("triqs_tprf::g_fk_t triqs_tprf::lattice_dyson_g0_fk (double mu, triqs_tprf::e_k_cvt e_k, triqs::mesh::refreq mesh, double delta)", doc = r"""Construct a non-interacting real frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(\omega, \mathbf{k})`

  Computes

  .. math::
     G^{(0)}_{a\bar{b}}(\omega, \mathbf{k}) = \left[
     (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
     \right]^{-1}_{a\bar{b}},

  using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, chemical potential :math:`\mu`,
  broadening :math:`\delta`, and a real frequency Green's function mesh.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

mesh
     real frequency mesh

delta
     broadening :math:`\delta`

Returns
-------
out
     Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_wk_cvt sigma_wk)", doc = r"""Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

 Computes

 .. math::
    G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
    (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a Matsubara frequency
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

module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_Dwk_cvt sigma_wk)")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::lattice_dyson_g_fk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_fk_cvt sigma_fk, double delta)", doc = r"""Construct an interacting real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`

 Computes

 .. math::
    G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega, \mathbf{k})
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, broadening :math:`\delta`, and a real frequency
 self energy :math:`\Sigma_{\bar{a}b}(\omega, \mathbf{k})`.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

sigma_fk
     real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega, \mathbf{k})`

delta
     broadening :math:`\delta`

Returns
-------
out
     real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::lattice_dyson_g_fk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_f_cvt sigma_f, double delta)", doc = r"""Construct an interacting real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`

 Computes

 .. math::
    G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
    (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, broadening :math:`\delta`, and a real frequency
 self energy :math:`\Sigma_{\bar{a}b}(\omega)`.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

sigma_fk
     real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega)`

delta
     broadening :math:`\delta`

Returns
-------
out
     real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`""")

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

module.add_function ("triqs_tprf::g_Dw_t triqs_tprf::lattice_dyson_g_w (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_Dw_cvt sigma_w)")

module.add_function ("triqs_tprf::g_f_t triqs_tprf::lattice_dyson_g_f (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_f_cvt sigma_f, double delta)", doc = r"""Construct an interacting real frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function :math:`G_{a\bar{b}}(\omega)`

 Computes

 .. math::
    G_{a\bar{b}}(\omega) = \frac{1}{N_k} \sum_\mathbf{k} \left[
    (i\omega_n + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
    \right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`,
 chemical potential :math:`\mu`, and a momentum independent real frequency
 self energy :math:`\Sigma_{\bar{a}b}(\omega)`.

Parameters
----------
mu
     chemical potential :math:`\mu`

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

sigma_f
     real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega)`

delta
     broadening :math:`\delta`

Returns
-------
out
     Real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wr_t triqs_tprf::fourier_wk_to_wr (triqs_tprf::g_wk_cvt g_wk)", doc = r"""Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(i\omega_n, \mathbf{k})\right\}`

Parameters
----------
g_wk
     k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

Returns
-------
out
     real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`""")

module.add_function ("triqs_tprf::g_Dwr_t triqs_tprf::fourier_wk_to_wr (triqs_tprf::g_Dwk_cvt g_wk)")

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

module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::fourier_wr_to_wk (triqs_tprf::g_Dwr_cvt g_wr)")

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

module.add_function ("triqs_tprf::g_Dtr_t triqs_tprf::fourier_wr_to_tr (triqs_tprf::g_Dwr_cvt g_wr, int nt = -1)")

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

module.add_function ("triqs_tprf::g_Dwr_t triqs_tprf::fourier_tr_to_wr (triqs_tprf::g_Dtr_cvt g_tr, int nw = -1)")

module.add_function ("triqs_tprf::g_fr_t triqs_tprf::fourier_fk_to_fr (triqs_tprf::g_fk_cvt g_fk)", doc = r"""Inverse fast fourier transform of real frequency Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(\omega, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(\omega, \mathbf{k})\right\}`

Parameters
----------
g_fk
     k-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`

Returns
-------
out
     real-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{r})`""")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::fourier_fr_to_fk (triqs_tprf::g_fr_cvt g_fr)", doc = r"""Fast fourier transform of real frequency Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(\omega, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(\omega, \mathbf{r}) \right\}`

Parameters
----------
g_fr
     real-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{r})`

Returns
-------
out
     k-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_Tr_t triqs_tprf::fourier_Tk_to_Tr (triqs_tprf::g_Tk_cvt g_Tk)", doc = r"""Inverse fast fourier transform of real time Green's function from k-space to real space

    Computes: :math:`G_{a\bar{b}}(t, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(t, \mathbf{k})\right\}`

Parameters
----------
g_Tk
     k-space real time Green's function :math:`G_{a\bar{b}}(T, \mathbf{k})`

Returns
-------
out
     real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`""")

module.add_function ("triqs_tprf::g_Tk_t triqs_tprf::fourier_Tr_to_Tk (triqs_tprf::g_Tr_cvt g_Tr)", doc = r"""Fast fourier transform of real time Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(t, \mathbf{r}) \right\}`

Parameters
----------
g_Tr
     real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`

Returns
-------
out
     k-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_Tk_t triqs_tprf::fourier_Tr_to_Tk (triqs_tprf::chi_Tr_cvt g_Tr)", doc = r"""Fast fourier transform of real time Green's function from real-space to k-space

    Computes: :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}c\bar{d}}(t, \mathbf{r}) \right\}`

Parameters
----------
g_Tr
     real-space real time Green's function :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{r})`

Returns
-------
out
     k-space real time Green's function :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::lindhard_chi00 (triqs_tprf::e_k_cvt e_k, triqs::mesh::imfreq mesh, double mu)", doc = r"""Generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`.

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
           U_{\bar{a}i}(\mathbf{k}) U^\dagger_{id}(\mathbf{k})
           U_{\bar{c}j}(\mathbf{k} + \mathbf{q}) U^\dagger_{jb}(\mathbf{k} + \mathbf{q})

    where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued
    dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

Parameters
----------
e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

mesh
     bosonic Matsubara frequency mesh

mu
     chemical potential :math:`\mu`

Returns
-------
out
     generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`

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


module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::lindhard_chi00 (triqs_tprf::e_k_cvt e_k, triqs::mesh::dlr_imfreq mesh, double mu)")

module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::lindhard_chi00 (triqs_tprf::e_k_cvt e_k, triqs::mesh::refreq mesh, double beta, double mu, double delta)", doc = r"""Generalized Lindhard susceptibility in the particle-hole channel and for real frequencies :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q})`.

    Analytic calculation of the generalized (non-interacting) Lindhard susceptibility
    in the particle-hole channel in real frequencies. The analytic expression is obtained using
    residue calculus to explicitly evaluate the matsubara sum of the fourier transformed imaginary
    time bubble product of two non-interacting single-particle Green's functions. The resulting
    expression is analytically continued to the real frequency axis.

    .. math::
       G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) =
       \left[ i\omega_n \cdot \mathbf{1} - \epsilon(\mathbf{k}) \right]^{-1} .

    Analytic continuation of the evaluation of the bubble diagram gives

    .. math::
        \chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q}) =
        \frac{1}{N_k} \sum_{\mathbf{k}} \sum_{ij}
        \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
            {\omega + i\delta + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
        \\ \times
        U_{\bar{a}i}(\mathbf{k}) U^\dagger_{id}(\mathbf{k})
        U_{\bar{c}j}(\mathbf{k} + \mathbf{q}) U^\dagger_{jb}(\mathbf{k} + \mathbf{q})

    where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued
    dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}

Parameters
----------
e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

mesh
     real frequency mesh

beta
     inverse temperature

mu
     chemical potential :math:`\mu`

delta
     broadening :math:`\delta`

Returns
-------
out
     real frequency generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::solve_rpa_PH (triqs_tprf::chi_wk_vt chi0, array_contiguous_view<std::complex<double>, 4> U)", doc = r"""Random Phase Approximation (RPA) in the particle-hole channel

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

module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::solve_rpa_PH (triqs_tprf::chi_fk_vt chi0, array_contiguous_view<std::complex<double>, 4> U)", doc = r"""Random Phase Approximation (RPA) in the particle-hole channel

     Computes the equation

     .. math::
         \chi(\bar{a}b\bar{c}d) = \big(
         \mathbb{1}
         - \chi^{(0)}(\bar{a}b\bar{B}A) U(A\bar{B}D\bar{C})
         \big)^{-1} \chi^{(0)}(\bar{C}D\bar{c}d)\,.

Parameters
----------
chi0
     bare particle-hole bubble :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{k}, \omega)`

U
     RPA static vertex as obtained from triqs_tprf.rpa_tensor.get_rpa_tensor :math:`U_{a\bar{b}c\bar{d}}`

Returns
-------
out
     RPA suceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, \omega)`""")

module.add_function ("g_w_t triqs_tprf::dlr_on_imfreq (triqs_tprf::g_Dc_cvt g_c, triqs::mesh::imfreq wmesh)")
module.add_function ("g_t_t triqs_tprf::dlr_on_imtime (triqs_tprf::g_Dc_cvt g_c, triqs::mesh::imtime tmesh)")
module.add_function ("chi_w_t triqs_tprf::dlr_on_imfreq (triqs_tprf::chi_Dc_cvt chi_c, triqs::mesh::imfreq wmesh)")
module.add_function ("chi_t_t triqs_tprf::dlr_on_imtime (triqs_tprf::chi_Dc_cvt chi_c, triqs::mesh::imtime tmesh)")

module.add_function ("std::tuple<chi_wk_t, chi_k_t> triqs_tprf::split_into_dynamic_wk_and_constant_k (triqs_tprf::chi_wk_cvt chi_wk)", doc = r"""Split mtrix-valued Green's function into dynamic and constant parts by tail fitting

Parameters
----------
chi_wk
     : general matrix-valued Green's function :math:`\chi(i\omega_n, \mathbf{k})`.

Returns
-------
out
     Tuple of chi_dyn_wk, the dynamic part of chi, which converges to zero for :math:`\omega_n \rightarrow \infty`, and chi_const_k, the part of chi that is constant in Matsubara frequency space :math:`\chi(\mathbf{k})`.""")


module.add_function ("triqs_tprf::g_fk_t triqs_tprf::add_dynamic_and_static(triqs_tprf::g_fk_t g_dyn_fk, triqs_tprf::e_k_t g_stat_k)", doc = r"""Add documentation!""")
module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::add_dynamic_and_static(triqs_tprf::chi_fk_t chi_dyn_fk, triqs_tprf::chi_k_t chi_stat_k)", doc = r"""Add documentation!""")
module.add_function ("triqs_tprf::g_wk_t triqs_tprf::add_dynamic_and_static(triqs_tprf::g_wk_t g_dyn_wk, triqs_tprf::e_k_t g_stat_k)", doc = r"""Add documentation!""")
module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::add_dynamic_and_static(triqs_tprf::chi_wk_t chi_dyn_wk, triqs_tprf::chi_k_t chi_stat_k)", doc = r"""Add documentation!""")
module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::add_dynamic_and_static(triqs_tprf::g_Dwk_t g_dyn_wk, triqs_tprf::e_k_t g_stat_k)", doc = r"""Add documentation!""")
module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::add_dynamic_and_static(triqs_tprf::chi_Dwk_t chi_dyn_wk, triqs_tprf::chi_k_t chi_stat_k)", doc = r"""Add documentation!""")

module.add_function ("double triqs_tprf::fermi(double e)", doc = r"""Add documentation!""")

module.add_function ("double triqs_tprf::bose(double e)", doc = r"""Add documentation!""")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::rho_k_from_g_wk (triqs_tprf::g_wk_cvt g_wk)", doc = r"""Density matrix from lattic Green's function
      
Parameters
----------
g_wk
     single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
rho_k
     density matrix :math:`\rho_{ab}(\mathbf{k})`
""")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::rho_k_from_g_wk (triqs_tprf::g_Dwk_cvt g_wk)")
                     
module.add_function ("triqs_tprf::g_wk_t triqs_tprf::gw_sigma (triqs_tprf::chi_wk_cvt W_wk, triqs_tprf::g_wk_cvt g_wk)", doc = r"""GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator for dynamic interactions

    Splits the interaction into a dynamic and a static part

    .. math ::
        W_{abcd}(i\omega_n, \mathbf{k}) =
            W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k})
            + V_{abcd}(\mathbf{k})

    by fitting the high-frequency tail.

    Fourier transforms the dynamic part of the interaction and the
    single-particle Green's function to imaginary time and real space.

    .. math::
        G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

    .. math::
        W^{(dyn)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
          \left\{ W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k}) \right\}

    computes the GW self-energy as the product

    .. math::
        \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) =
          - \sum_{cd} W^{(dyn)}_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

    and transforms back to frequency and momentum

    .. math::
        \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k}) =
          \mathcal{F} \left\{ \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) \right\}

    The self-energy of the static part of the interaction is calculated
    as the sum

    .. math::
        \Sigma^{(stat)}_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{\mathbf{q}} V_{abab}(\mathbf{k}) \rho(G(i\omega_n, \mathbf{k+q}))_{ab}

    where :math:`\rho(G(i\omega_n, \mathbf{k+q}))` is the density matrix of the
    single particle Green's function.

    The total GW self-energy is given by

    .. math::
        \Sigma_{ab}(i\omega_n, \mathbf{k}) =
          \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k})
          + \Sigma^{(stat)}_{ab}(\mathbf{k})

Parameters
----------
W_wk
     interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

g_wk
     single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
out
     GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::e_r_t triqs_tprf::hartree_sigma (triqs_tprf::chi_k_cvt v_k, triqs_tprf::e_r_cvt rho_r)")
module.add_function ("triqs_tprf::e_r_t triqs_tprf::fock_sigma (triqs_tprf::chi_r_cvt v_r, triqs_tprf::e_r_cvt rho_r)")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::hartree_sigma (triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_wk_cvt g_wk)", doc = r"""Hartree self energy :math:`\Sigma_{ab}(\mathbf{k})` calculator

    Computes the Hartree self-energy of a static interaction as the sum

    .. math::
        \Sigma_{ab}(\mathbf{k}) = \frac{1}{N_k}
          \sum_{\mathbf{q},cd} V_{abcd}(\mathbf{q}) \rho_{cd}(\mathbf{k} + \mathbf{q})

    where :math:`\rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k})` is the
    density matrix of the single particle Green's function.

Parameters
----------
V_k
     static interaction :math:`V_{abcd}(\mathbf{k})`

g_wk
     single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
out
     Hartree self-energy :math:`\Sigma_{ab}(\mathbf{k})`""")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::fock_sigma (triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_wk_cvt g_wk)", doc = r"""Fock self energy :math:`\Sigma_{ab}(\mathbf{k})` calculator 

    Computes the Fock self-energy of a static interaction as the sum

    .. math::
        \Sigma_{ab}(\mathbf{k}) = -\frac{1}{N_k}
          \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{q}) \rho_{dc}(\mathbf{k} + \mathbf{q})

    where :math:`\rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k})` is the density
    matrix of the single particle Green's function.

Parameters
----------
V_k
     static interaction :math:`V_{abcd}(\mathbf{k})`

g_wk
     single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
out
     Fock self-energy :math:`\Sigma_{ab}(\mathbf{k})`""")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::gw_sigma (triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_wk_cvt g_wk)", doc = r"""Static GW self energy :math:`\Sigma(\mathbf{k})` calculator

    Computes the GW self-energy (equivalent to the Fock self-energy)

Parameters
----------
V_k
     static interaction :math:`V_{abcd}(\mathbf{k})`

g_wk
     single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
out
     Static GW self-energy (Fock) :math:`\Sigma_{ab}(\mathbf{k})`""")

module.add_function ("triqs_tprf::g_tr_t triqs_tprf::gw_dynamic_sigma (triqs_tprf::chi_tr_cvt W_tr, triqs_tprf::g_tr_cvt g_tr)", doc = r"""Dynamic GW self energy :math:`\Sigma(\tau, \mathbf{r})` calculator

    Computes the GW self-energy as the product

    .. math::
        \Sigma_{ab}(\tau, \mathbf{r}) =
          - \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

Parameters
----------
W_tr
     interaction :math:`W_{abcd}(\tau, \mathbf{r})`

g_tr
     single particle Green's function :math:`G_{ab}(\tau, \mathbf{r})`

Returns
-------
out
     Dynamic GW self-energy :math:`\Sigma_{ab}(\tau, \mathbf{r})`""")

module.add_function ("triqs_tprf::g_Dtr_t triqs_tprf::gw_dynamic_sigma (triqs_tprf::chi_Dtr_cvt W_tr, triqs_tprf::g_Dtr_cvt g_tr)")
                     
module.add_function ("triqs_tprf::g_f_t triqs_tprf::g0w_dynamic_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta, mesh::brzone::value_t kpoint)", doc = r"""add documentation!""")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::g0w_dynamic_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta, mesh::brzone kmesh)", doc = r"""add documentation!""")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::g0w_dynamic_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta)", doc = r"""add documentation!""")

module.add_function ("array<std::complex<double>, 2> g0w_sigma(double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_k_cvt v_k, mesh::brzone::value_t kpoint)", doc = r"""Add some docs""")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::g0w_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_k_cvt v_k, mesh::brzone kmesh)", doc = r"""Add some docs""")

module.add_function ("triqs_tprf::e_k_t triqs_tprf::g0w_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_k_cvt v_k)", doc = r"""GW self energy :math:`\Sigma(\mathbf{k})` calculator for static interactions

    Computes the GW self-energy of a static interaction as the product

    .. math::
        \Sigma_{ab}(\mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
          U_{la}(\mathbf{k}+\mathbf{q}) U^\dagger_{bl}(\mathbf{k}+\mathbf{q})
          V_{abab}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l})

    where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued
    dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}

Parameters
----------
mu
     chemical potential :math:`\mu`

beta
     inverse temperature

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

V_k
     bare interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     static GW self-energy :math:`\Sigma_{ab}(\mathbf{k})`""")

module.add_function ("triqs_tprf::g_f_t triqs_tprf::g0w_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta, mesh::brzone::value_t kpoint)", doc = r"""add documentation!""")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::g0w_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta, mesh::brzone kmesh)", doc = r"""add documentation!""")

module.add_function ("triqs_tprf::g_fk_t triqs_tprf::g0w_sigma (double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta)", doc = r"""Real frequency GW self energy :math:`\Sigma(\omega, \mathbf{k})` calculator via the spectral representation

    Computes the spectral function of the dynamic part of the screened interaction

    .. math::
        W^{(spec)}_{ab}(\omega, \mathbf{k}) = \frac{-1}{\pi} \text{Im}
          \left( W_{abab}(\omega, \mathbf{k}) - V_{abab}(\mathbf{k}) \right)

    and constructs the GW self energy via the spectral representation

    .. math::
        \Sigma_{ab}(\omega, \mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
          U_{la}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{bl}(\mathbf{k}+\mathbf{q})
          V_{abab}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l}) \\
        + \frac{\delta_{\omega}}{N_k} \sum_{\mathbf{q}} \sum_{\omega'}
          U_{la}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{bl}(\mathbf{k}+\mathbf{q})
          W^{(spec)}_{ab}(\omega', \mathbf{q})
          \frac{n_B(\omega') + f(\epsilon_{\mathbf{k}+\mathbf{q}, l})}{\omega + i\delta + \omega' - \epsilon_{\mathbf{k}+\mathbf{q}, l} + \mu}

    where :math:`\delta_{\omega}` is the real-frequency mesh spacing and the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued
    dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

    .. math::
       \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
       = \delta_{ij} \epsilon_{\mathbf{k}, i}

Parameters
----------
mu
     chemical potential :math:`\mu`

beta
     inverse temperature

e_k
     discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

W_fk
     fully screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`

V_k
     bare interaction :math:`V_{abcd}(\mathbf{k})`

delta
     broadening :math:`\delta`

Returns
-------
out
     real frequency GW self-energy :math:`\Sigma_{ab}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::dynamical_screened_interaction_W (triqs_tprf::chi_wk_cvt PI_wk, triqs_tprf::chi_k_cvt V_k)", doc = r"""Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

Parameters
----------
PI_wk
     polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

V_k
     static bare interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::dynamical_screened_interaction_W (triqs_tprf::chi_Dwk_cvt PI_wk, triqs_tprf::chi_k_cvt V_k)")

module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::dynamical_screened_interaction_W (triqs_tprf::chi_fk_cvt PI_fk, triqs_tprf::chi_k_cvt V_k)", doc = r"""Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})`.

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) =
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \Pi_{fegh}(\omega, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(\omega, \mathbf{k})

Parameters
----------
PI_fk
     polarization bubble :math:`\Pi_{abcd}(\omega, \mathbf{k})`

V_k
     static bare interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::dynamical_screened_interaction_W (triqs_tprf::chi_wk_cvt PI_wk, triqs_tprf::chi_wk_cvt V_wk)", doc = r"""Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(i\omega_n, \mathbf{k})`.

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
        V_{abcd}(i\omega_n, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(i\omega_n, \mathbf{k}) \cdot
        \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

Parameters
----------
PI_wk
     polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`

V_wk
     dynamic bare interaction :math:`V_{abcd}(i\omega_n, \mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`""")


module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::dynamical_screened_interaction_W (triqs_tprf::chi_Dwk_cvt PI_wk, triqs_tprf::chi_Dwk_cvt V_wk)")

module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::dynamical_screened_interaction_W (triqs_tprf::chi_fk_cvt PI_fk, triqs_tprf::chi_fk_cvt V_fk)", doc = r"""Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(\omega, \mathbf{k})`.

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) =
        V_{abcd}(\omega, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(\omega, \mathbf{k}) \cdot
        \Pi_{fegh}(\omega, \mathbf{k}) \cdot
        W^{(full)}_{hgcd}(\omega, \mathbf{k})

Parameters
----------
PI_fk
     polarization bubble :math:`\Pi_{abcd}(\omega, \mathbf{k})`

V_fk
     dynamic bare interaction :math:`V_{abcd}(\omega, \mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility (triqs_tprf::chi_wk_cvt chi_wk, triqs_tprf::chi_k_cvt V_k)", doc = r"""Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})` and known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        V_{hgcd}(\mathbf{k})

Parameters
----------
chi_wk
     generalized susceptibility :math:`\chi_{abcd}(i\omega_n, \mathbf{k})`

V_k
     static bare interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility (triqs_tprf::chi_Dwk_cvt chi_wk, triqs_tprf::chi_k_cvt V_k)")

module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility (triqs_tprf::chi_fk_cvt chi_fk, triqs_tprf::chi_k_cvt V_k)", doc = r"""Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})` and known generalized susceptibility :math:`\chi(\omega, \mathbf{k})`

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) =
        V_{abcd}(\mathbf{k}) +
	    \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
        \chi_{fegh}(\omega, \mathbf{k}) \cdot
        V_{hgcd}(\mathbf{k})

Parameters
----------
chi_fk
     generalized susceptibility :math:`\chi_{abcd}(\omega, \mathbf{k})`

V_k
     static bare interaction :math:`V_{abcd}(\mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility (triqs_tprf::chi_wk_cvt chi_wk, triqs_tprf::chi_wk_cvt V_wk)", doc = r"""Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(i\omega_n, \mathbf{k})` and known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

    The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
        V_{abcd}(i\omega_n, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(i\omega_n, \mathbf{k}) \cdot
        \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
        V_{hgcd}(i\omega_n, \mathbf{k})

Parameters
----------
chi_wk
     generalized susceptibility :math:`\chi_{abcd}(i\omega_n, \mathbf{k})`

V_wk
     dynamic bare interaction :math:`V_{abcd}(i\omega_n, \mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`""")

module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility (triqs_tprf::chi_Dwk_cvt chi_wk, triqs_tprf::chi_Dwk_cvt V_wk)")

module.add_function ("triqs_tprf::chi_fk_t triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility (triqs_tprf::chi_fk_cvt chi_fk, triqs_tprf::chi_fk_cvt V_fk)", doc = r"""Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(\omega, \mathbf{k})` and known generalized susceptibility :math:`\chi(\omega, \mathbf{k})`

    The full screened interaction :math:`W(\omega, \mathbf{k})`
    is given by

    .. math::
        W^{(full)}_{abcd}(\omega, \mathbf{k}) =
        V_{abcd}(\omega, \mathbf{k}) +
	    \sum_{efgh} V_{abef}(\omega, \mathbf{k}) \cdot
        \chi_{fegh}(\omega, \mathbf{k}) \cdot
        V_{hgcd}(\omega, \mathbf{k})

Parameters
----------
chi_fk
     generalized susceptibility :math:`\chi_{abcd}(\omega, \mathbf{k})`

V_fk
     dynamic bare interaction :math:`V_{abcd}(\omega, \mathbf{k})`

Returns
-------
out
     dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::eliashberg_product (triqs_tprf::chi_wk_vt Gamma_pp, triqs_tprf::g_wk_vt g_wk, triqs_tprf::g_wk_vt delta_wk)", doc = r"""Linearized Eliashberg product via summation

     Computes the linearized Eliashberg product in the singlet/triplet channel given by

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k})
         =
         -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
         \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
         \\
         \times
         G_{c\bar{e}}(i\nu',\mathbf{k}')
         G_{d\bar{f}}(-i\nu',-\mathbf{k}')
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{f}}(i\nu',\mathbf{k}')\,,

     by summation.

Parameters
----------
Gamma_pp
     particle-particle vertex :math:`\Gamma^{\mathrm{s/t}}_{a\bar{b}c\bar{d}}(i\nu_n,\mathbf{k})`

g_wk
     single particle Green's function :math:`G_{a\bar{b}}(i\nu_n,\mathbf{k})`

delta_wk
     superconducting gap :math:`\Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})`

Returns
-------
out
     Gives the result of the product :math:`\Delta^{\mathrm{s/t}, \mathrm{out}}`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::eliashberg_product_fft (triqs_tprf::chi_tr_vt Gamma_pp_dyn_tr, triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_wk_vt g_wk, triqs_tprf::g_wk_vt delta_wk)", doc = r"""Linearized Eliashberg product via FFT

     Computes the linearized Eliashberg product in the singlet/triplet channel given by

     .. math::
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k})
        =
        -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
        \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
        \\
        \times
        G_{c\bar{e}}(i\nu',\mathbf{k}')
        G_{d\bar{f}}(-i\nu',-\mathbf{k}')
        \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{f}}(i\nu',\mathbf{k}')\,,

     by taking advantage of the convolution theorem.

     We therefore first calculate

     .. math::
         F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
         =
         G_{a\bar{c}}(i\nu,\mathbf{k})
         G_{b\bar{d}}(-i\nu,-\mathbf{k})
         \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})\,,

     which we then Fourier transform to imaginary time and real-space

     .. math::
        F^{\mathrm{s/t}}_{ab}(\tau,\mathbf{r})
        =
        \mathcal{F}^2
        \big(
        F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
        \big)\,.

     We then calculate first the dynamic gap

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
         =
         -\frac{1}{2}
         \Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})
         F^{\mathrm{s/t}}_{cd}(\tau, \mathbf{r})\,,

     and then the static gap

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})
         =
         -\frac{1}{2}
         \Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})
         F^{\mathrm{s/t}}_{cd}(\tau=0, \mathbf{r})\,.

     We then Fourier transform the dynamic gap to imaginary frequencies

     .. math::
         \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        =
        \mathcal{F}
        \big(
        \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
        \big)\,,

     and then add both component together

     .. math::
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        =
        \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        +
        \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})\,,

    and then finally Fourier transform to :math:`\mathbf{k}`-space

    .. math::
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})
        =
        \mathcal{F}
        \big(
        \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
        \big)\,.

Parameters
----------
Gamma_pp_dyn_tr
     dynamic part of the particle-particle vertex :math:`\Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})`

Gamma_pp_const_r
     static part of the particle-particle vertex :math:`\Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})`

g_wk
     one-particle Green's function :math:`G_{a\bar{b}}(i\nu_n,\mathbf{k})`

delta_wk
     superconducting gap :math:`\Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})`

Returns
-------
out
     Gives the result of the product :math:`\Delta^{\mathrm{s/t}, \mathrm{out}}`""")


module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::eliashberg_product_fft (triqs_tprf::chi_Dtr_vt Gamma_pp_dyn_tr, triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_Dwk_vt g_wk, triqs_tprf::g_Dwk_vt delta_wk)", doc = r"""Add documentation!""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::eliashberg_product_fft_constant (triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_wk_vt g_wk, triqs_tprf::g_wk_vt delta_wk)", doc = r"""""")

module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::eliashberg_product_fft_constant (triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_Dwk_vt g_wk, triqs_tprf::g_Dwk_vt delta_wk)", doc = r"""""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::eliashberg_g_delta_g_product (triqs_tprf::g_wk_vt g_wk, triqs_tprf::g_wk_vt delta_wk)", doc = r"""""")

module.add_function ("triqs_tprf::g_Dwk_t triqs_tprf::eliashberg_g_delta_g_product (triqs_tprf::g_Dwk_vt g_wk, triqs_tprf::g_Dwk_vt delta_wk)", doc = r"""""")

module.add_function ("std::tuple<chi_tr_t, chi_r_t> triqs_tprf::dynamic_and_constant_to_tr (triqs_tprf::chi_wk_vt Gamma_pp_dyn_wk, triqs_tprf::chi_k_vt Gamma_pp_const_k)", doc = r"""Fourier transform Gamma parts to imaginary time and real-space

Parameters
----------
Gamma_pp_dyn_wk
     : The dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`.

Gamma_pp_const_k
     : The part of Gamma that is constant in Matsubara frequency space :math:`\Gamma(\mathbf{k})`.

Returns
-------
out
     Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space, Gamma_pp_const_r, the constant part of Gamma in real-space.""")

module.add_function ("std::tuple<chi_Dtr_t, chi_r_t> triqs_tprf::dynamic_and_constant_to_tr (triqs_tprf::chi_Dwk_vt Gamma_pp_dyn_wk, triqs_tprf::chi_k_vt Gamma_pp_const_k)", doc = r"""Fourier transform Gamma parts to imaginary time and real-space

Parameters
----------
Gamma_pp_dyn_wk
     : The dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`.

Gamma_pp_const_k
     : The part of Gamma that is constant in Matsubara frequency space :math:`\Gamma(\mathbf{k})`.

Returns
-------
out
     Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space, Gamma_pp_const_r, the constant part of Gamma in real-space.""")


module.add_function ("triqs_tprf::e_r_t triqs_tprf::eliashberg_constant_gamma_f_product (triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_tr_t F_tr)", doc = r"""""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::construct_phi_wk (triqs_tprf::chi_wk_vt chi, array_contiguous_view<std::complex<double>, 4> U)", doc = r"""Computes reducible ladder vertex for the approximation of a local and static vertex.

    In this approximation the reducible ladder vertex in density/magnetic channel are given by

    .. math::
        \Phi^{\text{d/m}}_{a\overline{b}c\overline{d}}(Q)
        &\approx
        \frac{1}{(N_\mathbf{k}\beta)^2}
        \sum_{K'', K'''}
        U^{\text{d/m}}\chi^{\text{d/m}}(Q, K'', K''') U^{\text{d/m}}
        \\
        &\approx
        U^{\mathrm{d/m}}
        \chi^{\text{d/m}}(Q) U^{\mathrm{d/m}}\,,


    where all products are particle-hole products.
    The reducible ladder vertex in then only dependent on one bosonic frequency and momentum.
    It can then be used in :meth:`triqs_tprf.eliashberg.construct_gamma_singlet_rpa`
    or :meth:`triqs_tprf.eliashberg.construct_gamma__rpa` to construct the
    irreducible singlet/triplet vertex.

Parameters
----------
chi
     density/magnetic susceptibility  :math:`\chi^{\mathrm{d/m}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q})`

U
     density/magnetic local and static vertex  :math:`U^{\mathrm{d/m}}_{a\bar{b}c\bar{d}}`

Returns
-------
out
     The reducible ladder vertex in the density/magnetic channel :math:`\Phi^{\mathrm{d/m}}(i\omega_n,\mathbf{q})`""")

module.add_function ("array<std::complex<double>, 6> triqs_tprf::cluster_mesh_fourier_interpolation (array<double, 2> k_vecs, triqs_tprf::chi_wr_cvt chi)", doc = r"""""")

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

module.add_function ("triqs_tprf::chi_Dtr_t triqs_tprf::chi0_tr_from_grt_PH (triqs_tprf::g_Dtr_cvt g_tr)")
                     
module.add_function ("std::tuple<triqs_tprf::g_Tk_t, triqs_tprf::g_Tk_t> triqs_tprf::g0_Tk_les_gtr_from_e_k(triqs_tprf::e_k_cvt e_k, triqs::mesh::retime Tmesh, double beta)")

module.add_function ("triqs_tprf::chi_Tr_t triqs_tprf::chi0_Tr_from_g_Tr_PH (triqs_tprf::g_Tr_cvt g_Tr_les, triqs_tprf::g_Tr_cvt g_Tr_gtr)", doc = r"""Generalized susceptibility real time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r})`

  Computes

  .. math::
     \chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r}) =
     i G^<_{d\bar{a}}(t, \mathbf{r}) G^>_{b\bar{c}}(-t, -\mathbf{r})
     - i G^>_{d\bar{a}}(t, \mathbf{r}) G^<_{b\bar{c}}(-t, -\mathbf{r})

Parameters
----------
g_Tr_les
     Lesser real time Green's function in real-space, :math:`G^<_{a\bar{b}}(t, \mathbf{r})`.
g_Tr_gtr
     Greater real time Green's function in real-space, :math:`G^>_{a\bar{b}}(t, \mathbf{r})`.

Returns
-------
out
     Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r})` in real time and real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi0_wr_from_grt_PH (triqs_tprf::g_tr_cvt g_tr, int nw)", doc = r"""""")

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

module.add_function ("triqs_tprf::chi_Dwr_t triqs_tprf::chi_wr_from_chi_tr (triqs_tprf::chi_Dtr_cvt chi_tr, int nw)")

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

module.add_function ("triqs_tprf::chi_Dtr_t triqs_tprf::chi_tr_from_chi_wr (triqs_tprf::chi_Dwr_cvt chi_wr, int ntau = -1)")

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

module.add_function ("triqs_tprf::chi_Dwk_t triqs_tprf::chi_wk_from_chi_wr (triqs_tprf::chi_Dwr_cvt chi_wr)")

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

module.add_function ("triqs_tprf::chi_Dwr_t triqs_tprf::chi_wr_from_chi_wk (triqs_tprf::chi_Dwk_cvt chi_wk)")
                     
module.add_function ("chi_t_t::target_t::value_t triqs_tprf::chi_trapz_tau (triqs_tprf::chi_t_cvt chi_t)", doc = r"""""")

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

module.add_function ("triqs_tprf::chi_nr_t triqs_tprf::chi0_nr_from_gr_PH_at_specific_w (int nw_index, int nn, triqs_tprf::g_wr_cvt g_nr)", doc = r"""""")

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

module.add_function ("triqs_tprf::chi_kw_t triqs_tprf::chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH (triqs_tprf::chi_wnk_cvt chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn, triqs_tprf::chi_nn_cvt L_wn)", doc = r"""Dual lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

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

L_wn
     Local triangular particle-hole vertex function :math:`L^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu)`.

Returns
-------
out
     Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.""")

module.add_function ("gf<prod<triqs::mesh::brzone, triqs::mesh::imfreq>, tensor_valued<4>> triqs_tprf::chiq_sum_nu_from_g_wk_and_gamma_PH (triqs_tprf::gk_iw_t g_wk, triqs_tprf::g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1)", doc = r"""""")

module.add_function ("gf<prod<triqs::mesh::brzone, triqs::mesh::imfreq>, tensor_valued<4>> triqs_tprf::chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH (double mu, triqs_tprf::ek_vt e_k, triqs_tprf::g_iw_vt sigma_w, triqs_tprf::g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1)", doc = r"""""")

module.add_function ("gf<prod<triqs::mesh::brzone, triqs::mesh::imfreq>, tensor_valued<4>> triqs_tprf::chiq_sum_nu (triqs_tprf::chiq_t chiq)", doc = r"""""")

module.add_function ("gf<triqs::mesh::imfreq, tensor_valued<4>> triqs_tprf::chiq_sum_nu_q (triqs_tprf::chiq_t chiq)", doc = r"""""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::attatch_tri_vert (triqs_tprf::chi_nn_cvt L_wn, triqs_tprf::chi_kwnn_cvt chi_kwnn)", doc = r"""""")



module.generate_code()
