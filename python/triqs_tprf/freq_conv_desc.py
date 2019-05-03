# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/freq_conv.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m freq_conv -o freq_conv -C pytriqs --moduledoc="functionality for changing frequency conventions" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "freq_conv", doc = r"functionality for changing frequency conventions", app_name = "triqs_tprf")

# Imports
module.add_imports(*['pytriqs.gf'])

# Add here all includes
module.add_include("triqs_tprf/freq_conv.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", doc = r"""Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g_iw_t triqs_tprf::block_iw_AB_to_matrix_valued (triqs_tprf::b_g_iw_vt bg_AB)", doc = r"""""")

module.add_function ("void triqs_tprf::block_3nu_AABB_to_tensor_valued (triqs_tprf::b_g2_iw_vt bg2_AABB, triqs_tprf::g2_iw_vt g2)", doc = r"""""")

module.add_function ("void triqs_tprf::get_magnetic_component (triqs_tprf::g2_iw_vt g2, triqs_tprf::g2_iw_vt g2_m)", doc = r"""""")

module.add_function ("void triqs_tprf::from_3nu_PH (triqs_tprf::g2_iw_vt g2_ch, triqs_tprf::g2_iw_vt g2)", doc = r"""""")

module.add_function ("void triqs_tprf::from_3nu_PH_bar (triqs_tprf::g2_iw_vt g2_ch, triqs_tprf::g2_iw_vt g2)", doc = r"""""")

module.add_function ("void triqs_tprf::from_3nu_PP (triqs_tprf::g2_iw_vt g2_ch, triqs_tprf::g2_iw_vt g2)", doc = r"""""")



module.generate_code()