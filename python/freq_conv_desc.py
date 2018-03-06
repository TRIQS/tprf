
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "freq_conv", doc = "Two-particle Green's function freqiency notation converters", app_name = "freq_conv")

# All the triqs C++/Python modules
import pytriqs.gf

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/freq_conv.hpp") # Manually added

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
#include <cpp2py/converters/variant.hpp>

using namespace triqs::gfs;
using namespace triqs::lattice;
using namespace tprf;
""")

module.add_function ("gf<imfreq, matrix_valued> block_iw_AB_to_matrix_valued(block_gf_view<imfreq, matrix_valued> bg_AB)", doc = """""")

module.add_function ("void block_3nu_AABB_to_tensor_valued(block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> bg2_AABB, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")

module.add_function ("void get_magnetic_component(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2,  gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_m)", doc = """""")

module.add_function ("void from_3nu_PH(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_ch, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")
module.add_function ("void from_3nu_PH_bar(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_ch, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")
module.add_function ("void from_3nu_PP(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_ch, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")

module.generate_code()
