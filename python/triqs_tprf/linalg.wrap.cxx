
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// ==================== module functions ====================

// identity_PH
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt g) { return triqs_tprf::identity_PH(g); }, "g"),
                                                        c2py::cfun([](triqs_tprf::g2_nn_vt g) { return triqs_tprf::identity_PH(g); }, "g")};

// identity_PH_bar
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt g) { return triqs_tprf::identity_PH_bar(g); }, "g"),
                                                        c2py::cfun([](triqs_tprf::g2_nn_vt g) { return triqs_tprf::identity_PH_bar(g); }, "g")};

// identity_PP
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt g) { return triqs_tprf::identity_PP(g); }, "g"),
                                                        c2py::cfun([](triqs_tprf::g2_nn_vt g) { return triqs_tprf::identity_PP(g); }, "g")};

// inverse_PH
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt g) { return triqs_tprf::inverse_PH(g); }, "g"),
                                                        c2py::cfun([](triqs_tprf::g2_nn_vt g) { return triqs_tprf::inverse_PH(g); }, "g")};

// inverse_PH_bar
static auto const _c2py_fun_4 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt g) { return triqs_tprf::inverse_PH_bar(g); }, "g"),
                                                        c2py::cfun([](triqs_tprf::g2_nn_vt g) { return triqs_tprf::inverse_PH_bar(g); }, "g")};

// inverse_PP
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt g) { return triqs_tprf::inverse_PP(g); }, "g"),
                                                        c2py::cfun([](triqs_tprf::g2_nn_vt g) { return triqs_tprf::inverse_PP(g); }, "g")};

// product_PH
static auto const _c2py_fun_6 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B) { return triqs_tprf::product_PH(A, B); }, "A", "B"),
                           c2py::cfun([](triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B) { return triqs_tprf::product_PH(A, B); }, "A", "B")};

// product_PH_bar
static auto const _c2py_fun_7 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B) { return triqs_tprf::product_PH_bar(A, B); }, "A", "B"),
                           c2py::cfun([](triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B) { return triqs_tprf::product_PH_bar(A, B); }, "A", "B")};

// product_PP
static auto const _c2py_fun_8 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B) { return triqs_tprf::product_PP(A, B); }, "A", "B"),
                           c2py::cfun([](triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B) { return triqs_tprf::product_PP(A, B); }, "A", "B")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
Two-particle response-function identity operator :math:`\mathbf{1}` in the particle-hole channel (PH).

Constructs the unity-operator in the given channel

.. math::

   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

Storage is allocated and the result is returned by value.

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
g : {par_0}
   two-particle response function :math:`g \equiv g_{abcd}(\omega, \nu, \nu')` determining the shape and size of the unity operator

Returns
-------
{ret_0}
   the unity operator :math:`\mathbf{1}`, in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}}, {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC(
Two-particle response-function identity operator :math:`\mathbf{1}` in the particle-hole-bar channel (PH-bar).

Constructs the unity-operator in the given channel

.. math::

   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

Storage is allocated and the result is returned by value.

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
g : {par_0}
   two-particle response function :math:`g \equiv g_{abcd}(\omega, \nu, \nu')` determining the shape and size of the unity operator

Returns
-------
{ret_0}
   the unity operator :math:`\mathbf{1}`, in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}}, {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC(
Two-particle response-function identity operator :math:`\mathbf{1}` in the particle-particle channel (PP).

Constructs the unity-operator in the given channel

.. math::

   \mathbf{1}_{abcd}(\omega,\nu,\nu') =
   \mathbf{1}_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega)
   \equiv
   \delta_{\nu\nu'} \delta_{\alpha\gamma} \delta_{\beta\delta}

where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

Storage is allocated and the result is returned by value.

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
g : {par_0}
   two-particle response function :math:`g \equiv g_{abcd}(\omega, \nu, \nu')` determining the shape and size of the unity operator

Returns
-------
{ret_0}
   the unity operator :math:`\mathbf{1}`, in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}}, {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
Two-particle response-function inversion :math:`[g]^{-1}` in the particle-hole channel (PH).

The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')` 
is cast to matrix form and inverted

.. math::

   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

Storage is allocated and the inverse is returned by value.

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
g : {par_0}
   two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   :math:`[g]^{-1}` in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}}, {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC(
Two-particle response-function inversion :math:`[g]^{-1}` in the particle-hole-bar channel (PH-bar).

The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')` 
is cast to matrix form and inverted

.. math::

   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

Storage is allocated and the inverse is returned by value.

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
g : {par_0}
   two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   :math:`[g]^{-1}` in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}}, {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_5 = _c2py_fun_5.doc(R"DOC(
Two-particle response-function inversion :math:`[g]^{-1}` in the particle-particle channel (PP).

The two-particle response function :math:`g_{abcd}(\omega, \nu, \nu')` 
is cast to matrix form and inverted

.. math::

   [g]^{-1} = [ g_{\{\nu\alpha\beta\}, \{\nu'\gamma\delta\}}(\omega) ]^{-1}

where the mapping of target-space indices :math:`\{a, b, c, d \}` to :math:`\{\alpha, \beta\}, \{\gamma, \delta\}` is channel dependent.

Storage is allocated and the inverse is returned by value.

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
g : {par_0}
   two-particle response function to invert, :math:`g \equiv g_{abcd}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   :math:`[g]^{-1}` in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}}, {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_6 = _c2py_fun_6.doc(R"DOC(
Two-particle response-function product :math:`A * B` in the particle-hole channel (PH).

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

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
A : {par_0}
   two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
B : {par_1}
   two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   :math:`(A * B)` in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_7 = _c2py_fun_7.doc(R"DOC(
Two-particle response-function product :math:`A * B` in the particle-hole-bar channel (PH-bar).

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

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
A : {par_0}
   two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
B : {par_1}
   two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   :math:`(A * B)` in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_8 = _c2py_fun_8.doc(R"DOC(
Two-particle response-function product :math:`A * B` in the particle-particle channel (PP).

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

.. note::

   Assign to gf (g2_iw_t) yields move operation while assigning to gf_view (g2_iw_vt) causes extra copy operation

Parameters
----------
A : {par_0}
   two-particle response function :math:`A \equiv A_{abcd}(\omega, \nu, \nu')`
B : {par_1}
   two-particle response function :math:`B \equiv A_{abcd}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   :math:`(A * B)` in the given channel
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g2_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"identity_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"identity_PH_bar", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"identity_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"inverse_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"inverse_PH_bar", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"inverse_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"product_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {"product_PH_bar", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"product_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "linalg",                                                                           /* name of module */
   R"RAWDOC(Product, Inverse and Identity for two-particle response functions)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_linalg() {

  if (not c2py::check_python_version("linalg")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)

#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
