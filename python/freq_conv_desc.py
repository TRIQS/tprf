# Generated automatically using the command :
# c++2py.py -mpytriqs.applications.tprf.chi_from_gg2 ../c++/chi_from_gg2.hpp -I../c++ --moduledoc "Green's function generalized susceptibility tools"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.tprf.freq_conv", doc = "Green's function generalized susceptibility tools", app_name = "pytriqs.applications.tprf.freq_conv")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/freq_conv.hpp") # Manually added

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

module.add_function ("gf<imfreq, matrix_valued> block_iw_AB_to_matrix_valued(block_gf_view<imfreq, matrix_valued> bg_AB)", doc = """""")

module.add_function ("void block_3nu_AABB_to_tensor_valued(block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> bg2_AABB, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")

module.add_function ("void get_magnetic_component(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2,  gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_m)", doc = """""")

module.add_function ("void from_3nu_PH(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_ch, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")
module.add_function ("void from_3nu_PH_bar(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_ch, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")
module.add_function ("void from_3nu_PP(gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2_ch, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> g2)", doc = """""")

module.generate_code()
