#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_lattice_GUARDS
#define C2PY_HXX_DECLARATION_lattice_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_tprf::Channel_t> = true;
template <>
const std::map<triqs_tprf::Channel_t, str_t> c2py::enum_to_string<triqs_tprf::Channel_t> = {{triqs_tprf::Channel_t::PP, "PP"},
                                                                                            {triqs_tprf::Channel_t::PH, "PH"},
                                                                                            {triqs_tprf::Channel_t::PH_bar, "PH_bar"}};
#endif