
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "chi_from_gg2", doc = "Green's function generalized susceptibility tools", app_name = "chi_from_gg2")

# All the triqs C++/Python modules
import pytriqs.gf
import pytriqs.lattice

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/chi_from_gg2.hpp") # Manually added

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/lattice/brillouin_zone.hpp>

#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <cpp2py/converters/variant.hpp>

using namespace triqs::gfs;
using namespace triqs::lattice;
using namespace tprf;
""")

module.add_function ("g2_iw_t chi0_from_gg2_PH(gf_view<imfreq, matrix_valued> g, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")
module.add_function ("g2_iw_t chi0_from_gg2_PP(gf_view<imfreq, matrix_valued> g, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")

module.add_function ("g2_iw_t chi_from_gg2_PH(gf_view<imfreq, matrix_valued> g, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")
module.add_function ("g2_iw_t chi_from_gg2_PP(gf_view<imfreq, matrix_valued> g, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")

module.generate_code()
