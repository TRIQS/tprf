
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

// chi_ph_magnetic
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](int nw, int nwf, double beta, double U) { return triqs_tprf::hubbard_atom::chi_ph_magnetic(nw, nwf, beta, U); }, "nw", "nwf", "beta", "U")};

// gamma_ph_magnetic
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](int nw, int nwf, double beta, double U) { return triqs_tprf::hubbard_atom::gamma_ph_magnetic(nw, nwf, beta, U); }, "nw", "nwf", "beta", "U")};

// single_particle_greens_function
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](int nw, double beta, double U) { return triqs_tprf::hubbard_atom::single_particle_greens_function(nw, beta, U); }, "nw", "beta", "U")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(
   R"DOC(
Magnetic susceptibility of the Hubbard atom at half-filling :math:`\chi(\omega, \nu, \nu')`

    Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
    please cite the paper if you use this function!

Parameters
----------
nw : {par_0}
   number of bosonic Matsubara frequencies
nwf : {par_1}
   number of fermionic Matsubara frequencies
beta : {par_2}
   inverse temperature
U : {par_3}
   Hubbard U interaction parmeter

Returns
-------
{ret_0}
   chi magnetic susceptibility :math:`\chi(\omega, \nu, \nu')`
)DOC",
   {{c2py::python_typename<int>()}, {c2py::python_typename<int>()}, {c2py::python_typename<double>()}, {c2py::python_typename<double>()}},
   {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(
   R"DOC(
Magnetic vertex function in the particle-hole channel of the Hubbard atom at half-filling :math:`\Gamma(\omega, \nu, \nu')`

   Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
   please cite the paper if you use this function!

Parameters
----------
nw : {par_0}
   number of bosonic Matsubara frequencies
nwf : {par_1}
   number of fermionic Matsubara frequencies
beta : {par_2}
   inverse temperature
U : {par_3}
   Hubbard U interaction parmeter

Returns
-------
{ret_0}
   gamma magnetic susceptibility :math:`\Gamma(\omega, \nu, \nu')`
)DOC",
   {{c2py::python_typename<int>()}, {c2py::python_typename<int>()}, {c2py::python_typename<double>()}, {c2py::python_typename<double>()}},
   {c2py::python_typename<triqs_tprf::g2_iw_t>()});
static const auto _c2py_doc_2 =
   _c2py_fun_2.doc(R"DOC(
Single-particle Green's function of the Hubbard atom at half-filling

    Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
    please cite the paper if you use this function!

    

.. math::

   G(i\omega_n) = \frac{1}{i\omega_n - \frac{U^2}{4 i\omega_n}}

Parameters
----------
nw : {par_0}
   number of Matsubara frequencies
beta : {par_1}
   inverse temperature
U : {par_2}
   Hubbard U interaction parmeter

Returns
-------
{ret_0}
   g single-particle Green's function of the Hubbard atom :math:`G(i\omega_n)`
)DOC",
                   {{c2py::python_typename<int>()}, {c2py::python_typename<double>()}, {c2py::python_typename<double>()}},
                   {c2py::python_typename<triqs_tprf::g_iw_t>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"chi_ph_magnetic", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"gamma_ph_magnetic", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"single_particle_greens_function", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "hubbard_atom",                                                    /* name of module */
                                        R"RAWDOC(Exact correlation functions for the hubbard atom)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_hubbard_atom() {

  if (not c2py::check_python_version("hubbard_atom")) return NULL;

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
