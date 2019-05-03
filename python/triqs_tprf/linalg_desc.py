# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/linalg.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m linalg -o linalg -C pytriqs --moduledoc="Product, Inverse and Identity for two-particle response functions" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "linalg", doc = r"Product, Inverse and Identity for two-particle response functions", app_name = "triqs_tprf")

# Imports
module.add_imports(*['pytriqs.gf'])

# Add here all includes
module.add_include("triqs_tprf/linalg.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", doc = r"""Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::inverse_PH (triqs_tprf::g2_iw_vt g)", doc = r"""Two-particle response-function inversion :math:`[g]^{-1}` in the particle-hole channel (PH).

 The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')`
 is cast to matrix form and inverted

 .. math::
   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the inverse is returned by value.

Parameters
----------
g
     two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`

Returns
-------
out
     :math:`[g]^{-1}` in the given channel""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::inverse_PP (triqs_tprf::g2_iw_vt g)", doc = r"""Two-particle response-function inversion :math:`[g]^{-1}` in the particle-particle channel (PP).

 The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')`
 is cast to matrix form and inverted

 .. math::
   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the inverse is returned by value.

Parameters
----------
g
     two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`

Returns
-------
out
     :math:`[g]^{-1}` in the given channel""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::inverse_PH_bar (triqs_tprf::g2_iw_vt g)", doc = r"""Two-particle response-function inversion :math:`[g]^{-1}` in the particle-hole-bar channel (PH-bar).

 The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')`
 is cast to matrix form and inverted

 .. math::
   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the inverse is returned by value.

Parameters
----------
g
     two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`

Returns
-------
out
     :math:`[g]^{-1}` in the given channel""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::inverse_PH (triqs_tprf::g2_nn_vt g)", doc = r"""""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::inverse_PP (triqs_tprf::g2_nn_vt g)", doc = r"""""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::inverse_PH_bar (triqs_tprf::g2_nn_vt g)", doc = r"""""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::product_PH (triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B)", doc = r"""Two-particle response-function product :math:`A * B` in the particle-hole channel (PH).

 The two-particle response functions :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
 and :math:`B \equiv B_{abcd}(\omega, \nu, \nu')` are cast to matrix form and their
 product is computed

 .. math::
   (A * B)_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   = \sum_{\bar{\nu}ab}
   A_{\{\nu\alpha\beta\}, \{\bar{\nu}ab\}}(\omega)
   B_{\{\bar{\nu}ab\}, \{\nu'\gamma\delta\}}(\omega)

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the product is returned by value.

Parameters
----------
A
     two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`

B
     two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`

Returns
-------
out
     :math:`(A * B)` in the given channel""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::product_PP (triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B)", doc = r"""Two-particle response-function product :math:`A * B` in the particle-particle channel (PP).

 The two-particle response functions :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
 and :math:`B \equiv B_{abcd}(\omega, \nu, \nu')` are cast to matrix form and their
 product is computed

 .. math::
   (A * B)_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   = \sum_{\bar{\nu}ab}
   A_{\{\nu\alpha\beta\}, \{\bar{\nu}ab\}}(\omega)
   B_{\{\bar{\nu}ab\}, \{\nu'\gamma\delta\}}(\omega)

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the product is returned by value.

Parameters
----------
A
     two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`

B
     two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`

Returns
-------
out
     :math:`(A * B)` in the given channel""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::product_PH_bar (triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B)", doc = r"""Two-particle response-function product :math:`A * B` in the particle-hole-bar channel (PH-bar).

 The two-particle response functions :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
 and :math:`B \equiv B_{abcd}(\omega, \nu, \nu')` are cast to matrix form and their
 product is computed

 .. math::
   (A * B)_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   = \sum_{\bar{\nu}ab}
   A_{\{\nu\alpha\beta\}, \{\bar{\nu}ab\}}(\omega)
   B_{\{\bar{\nu}ab\}, \{\nu'\gamma\delta\}}(\omega)

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the product is returned by value.

Parameters
----------
A
     two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`

B
     two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`

Returns
-------
out
     :math:`(A * B)` in the given channel""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::product_PH (triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B)", doc = r"""""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::product_PP (triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B)", doc = r"""""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::product_PH_bar (triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B)", doc = r"""""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::identity_PH (triqs_tprf::g2_iw_vt g)", doc = r"""Two-particle response-function identity operator :math:`\mathbf{1}` in the particle-hole channel (PH).

 Constructs the unity-operator in the given channel

 .. math::
   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the result is returned by value.

Parameters
----------
A
     two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')` determinig the shape and size of the unity operator

Returns
-------
out
     the unity operator :math:`\mathbf{1}`, in the given channel""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::identity_PP (triqs_tprf::g2_iw_vt g)", doc = r"""Two-particle response-function identity operator :math:`\mathbf{1}` in the particle-particle channel (PP).

 Constructs the unity-operator in the given channel

 .. math::
   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the result is returned by value.

Parameters
----------
A
     two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')` determinig the shape and size of the unity operator

Returns
-------
out
     the unity operator :math:`\mathbf{1}`, in the given channel""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::identity_PH_bar (triqs_tprf::g2_iw_vt g)", doc = r"""Two-particle response-function identity operator :math:`\mathbf{1}` in the particle-hole-bar channel (PH-bar).

 Constructs the unity-operator in the given channel

 .. math::
   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

 where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

 Storage is allocated and the result is returned by value.

Parameters
----------
A
     two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')` determinig the shape and size of the unity operator

Returns
-------
out
     the unity operator :math:`\mathbf{1}`, in the given channel""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::identity_PH (triqs_tprf::g2_nn_vt g)", doc = r"""""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::identity_PP (triqs_tprf::g2_nn_vt g)", doc = r"""""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::identity_PH_bar (triqs_tprf::g2_nn_vt g)", doc = r"""""")



module.generate_code()