
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

// chi0_from_gg2_PH
static auto const _c2py_fun_0 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::chi0_from_gg2_PH(g, g2); }, "g", "g2")};

// chi0_from_gg2_PP
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::chi0_from_gg2_PP(g, g2); }, "g", "g2")};

// chi0_tau_from_g_tau_PH
static auto const _c2py_fun_2 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_tau_cvt g) { return triqs_tprf::chi0_tau_from_g_tau_PH(g); }, "g")};

// chi_from_gg2_PH
static auto const _c2py_fun_3 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::chi_from_gg2_PH(g, g2); }, "g", "g2")};

// chi_from_gg2_PP
static auto const _c2py_fun_4 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::chi_from_gg2_PP(g, g2); }, "g", "g2")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Hole channel

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
   - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega + \nu)

Parameters
----------
g : {par_0}
   single particle Green's function :math:`G_{ab}(\nu)`
g2 : {par_1}
   two-particle Green's function with the mesh to use for
     :math:`\chi^{(0)}`

Returns
-------
{ret_0}
   chi0 particle-hole bubble
     :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC(
Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Particle channel

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
   - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega - \nu)

Parameters
----------
g : {par_0}
   single particle Green's function :math:`G_{ab}(\nu)`
g2 : {par_1}
   two-particle Green's function with the mesh to use for
     :math:`\chi^{(0)}`

Returns
-------
{ret_0}
   chi0 particle-particle bubble
     :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_2 =
   _c2py_fun_2.doc(R"DOC(
Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Hole channel

    Computes

    

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau) =
   G_{da}(\nu) \cdot G_{bc}(\beta - \tay)

Parameters
----------
g : {par_0}
   single particle Green's function :math:`G_{ab}(\tau)`

Returns
-------
{ret_0}
   chi0 particle-hole bubble
     :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau)`
)DOC",
                   {{c2py::python_typename<triqs_tprf::g_tau_cvt>()}}, {c2py::python_typename<triqs_tprf::chi2_tau_t>()});
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
Generalized susceptibility :math:`\chi = G^{(2)} - GG` in the
Particle-Hole channel

    Computes

    

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
   G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   - \beta \delta_{\omega} G_{ba}(\nu) \cdot G_{dc}(\nu')

Parameters
----------
g : {par_0}
   single particle Green's function :math:`G_{ab}(\nu)`
g2 : {par_1}
   two-particle Green's function
     :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   chi generalized particle-hole susceptibility
     :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC(
Generalized susceptibility :math:`\chi = G^{(2)} - GG` in the
Particle-Particle channel

    Computes

    

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
   G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
   - \beta \delta_{\nu + \nu' - \omega} G_{ba}(\nu) \cdot G_{dc}(\nu')

Parameters
----------
g : {par_0}
   single particle Green's function :math:`G_{ab}(\nu)`
g2 : {par_1}
   two-particle Green's function
     :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`

Returns
-------
{ret_0}
   chi generalized particle-hole susceptibility
     :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_iw_vt>()}, {c2py::python_typename<triqs_tprf::g2_iw_vt>()}},
                                                {c2py::python_typename<triqs_tprf::g2_iw_t>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"chi0_from_gg2_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"chi0_from_gg2_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"chi0_tau_from_g_tau_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"chi_from_gg2_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"chi_from_gg2_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "chi_from_gg2",                                               /* name of module */
                                        R"RAWDOC(Calculation of generalized susceptibilities)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_chi_from_gg2() {

  if (not c2py::check_python_version("chi_from_gg2")) return NULL;

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
