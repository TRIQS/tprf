
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "linalg", doc = "Green's function linear algebra tools", app_name = "linalg")

# All the triqs C++/Python modules
import pytriqs.gf

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/linalg.hpp") # Manually added

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/variant.hpp>

using namespace triqs::gfs;
using namespace triqs::lattice;
using namespace tprf;
""")

module.add_function ("g2_iw_t identity_PH(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g)", doc = """""")
module.add_function ("g2_iw_t identity_PP(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g)", doc = """""")
module.add_function ("g2_iw_t identity_PH_bar(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g)", doc = """""")

module.add_function ("g2_iw_t product_PH(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> A, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> B)", doc = """""")
module.add_function ("g2_iw_t product_PP(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> A, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> B)", doc = """""")
module.add_function ("g2_iw_t product_PH_bar(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> A, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> B)", doc = """""")

module.add_function ("g2_iw_t inverse_PH(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g)", doc = """""")
module.add_function ("g2_iw_t inverse_PP(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g)", doc = """""")
module.add_function ("g2_iw_t inverse_PH_bar(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g)", doc = """""")

module.generate_code()
