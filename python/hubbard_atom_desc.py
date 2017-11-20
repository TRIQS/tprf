
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "hubbard_atom", doc = "Hubbard atom response functions", app_name = "hubbard_atom")

# All the triqs C++/Python modules
import pytriqs.gf

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/hubbard_atom.hpp") # Manually added

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
using namespace tprf::hubbard_atom;
using namespace tprf;
""")

module.add_function("g_iw_t single_particle_greens_function(int nw, double beta, double U)", doc = """""")
module.add_function("g2_iw_t chi_ph_magnetic(int nw, int nwf, double beta, double U)", doc = """""")

module.generate_code()
