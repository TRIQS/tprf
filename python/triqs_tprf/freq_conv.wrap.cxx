
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

// block_3nu_AABB_to_tensor_valued
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::b_g2_iw_vt bg2_AABB, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::block_3nu_AABB_to_tensor_valued(bg2_AABB, g2); },
              "bg2_AABB", "g2")};

// block_iw_AB_to_matrix_valued
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::b_g_iw_vt bg_AB) { return triqs_tprf::block_iw_AB_to_matrix_valued(bg_AB); }, "bg_AB")};

// from_3nu_PH
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g2_iw_vt g2_ch, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::from_3nu_PH(g2_ch, g2); }, "g2_ch", "g2")};

// from_3nu_PH_bar
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g2_iw_vt g2_ch, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::from_3nu_PH_bar(g2_ch, g2); }, "g2_ch", "g2")};

// from_3nu_PP
static auto const _c2py_fun_4 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g2_iw_vt g2_ch, triqs_tprf::g2_iw_vt g2) { return triqs_tprf::from_3nu_PP(g2_ch, g2); }, "g2_ch", "g2")};

// get_magnetic_component
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g2_iw_vt g2, triqs_tprf::g2_iw_vt g2_m) { return triqs_tprf::get_magnetic_component(g2, g2_m); }, "g2", "g2_m")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC()DOC");
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC()DOC");
static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC()DOC");
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC()DOC");
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC()DOC");
static const auto _c2py_doc_5 = _c2py_fun_5.doc(R"DOC()DOC");
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"block_3nu_AABB_to_tensor_valued", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"block_iw_AB_to_matrix_valued", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"from_3nu_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"from_3nu_PH_bar", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"from_3nu_PP", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"get_magnetic_component", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "freq_conv",                                                       /* name of module */
                                        R"RAWDOC(functionality for changing frequency conventions)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_freq_conv() {

  if (not c2py::check_python_version("freq_conv")) return NULL;

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
