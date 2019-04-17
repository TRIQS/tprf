# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/linalg.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m linalg -o linalg -C pytriqs --moduledoc="Product, Inverse and Identity for two-particle response functions" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "linalg", doc = "Product, Inverse and Identity for two-particle response functions", app_name = "triqs_tprf")

# Imports
import pytriqs.gf

# Add here all includes
module.add_include("triqs_tprf/linalg.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", """Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::inverse_PH (triqs_tprf::g2_iw_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::inverse_PP (triqs_tprf::g2_iw_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::inverse_PH_bar (triqs_tprf::g2_iw_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::inverse_PH (triqs_tprf::g2_nn_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::inverse_PP (triqs_tprf::g2_nn_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::inverse_PH_bar (triqs_tprf::g2_nn_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::product_PH (triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::product_PP (triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::product_PH_bar (triqs_tprf::g2_iw_vt A, triqs_tprf::g2_iw_vt B)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::product_PH (triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::product_PP (triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::product_PH_bar (triqs_tprf::g2_nn_vt A, triqs_tprf::g2_nn_vt B)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::identity_PH (triqs_tprf::g2_iw_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::identity_PP (triqs_tprf::g2_iw_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::identity_PH_bar (triqs_tprf::g2_iw_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::identity_PH (triqs_tprf::g2_nn_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::identity_PP (triqs_tprf::g2_nn_vt g)", doc = """""")

module.add_function ("triqs_tprf::g2_nn_t triqs_tprf::identity_PH_bar (triqs_tprf::g2_nn_vt g)", doc = """""")



module.generate_code()