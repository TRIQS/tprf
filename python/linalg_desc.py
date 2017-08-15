# Generated automatically using the command :
# c++2py.py -mpytriqs.applications.tprf.linalg ../c++/linalg.hpp -I../c++ --moduledoc "Green's function linear algebra tools"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.tprf.linalg", doc = "Green's function linear algebra tools", app_name = "pytriqs.applications.tprf.linalg")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/linalg.hpp") # Manually added

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/optional.hpp>
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/set.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/optional.hpp>
#include <triqs/python_tools/converters/variant.hpp>
#include <triqs/python_tools/converters/tuple.hpp>
#include <triqs/python_tools/converters/arrays.hpp>
#include <triqs/python_tools/converters/gf.hpp>
#include <triqs/python_tools/converters/block_gf.hpp>
#include <triqs/python_tools/converters/block2_gf.hpp>
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
