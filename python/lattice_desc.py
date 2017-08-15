# Generated automatically using the command :
# c++2py.py -mpytriqs.applications.tprf.lattice ../c++/lattice.hpp -I../c++ --moduledoc "Green's function lattice tools"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.tprf.lattie", doc = "Green's function lattice tools", app_name = "pytriqs.applications.tprf.lattice")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/types.hpp") # Manually added
module.add_include("../c++/lattice.hpp") # Manually added

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/lattice/brillouin_zone.hpp>
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

#module.add_function ("gk_iw_t g0k_from_ek(double mu, ek_cvt ek, g_iw_t::mesh_t mesh)", doc = """""")
#module.add_function ("gk_iw_t gk_from_ek_sigma(double mu, ek_cvt ek, g_iw_cvt sigma)", doc = """""")
module.add_function ("gr_iw_t gr_from_gk(gk_iw_vt gk)", doc = """""")

module.generate_code()
