
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "lattice", doc = "Green's function lattice tools", app_name = "lattice")

# All the triqs C++/Python modules
import pytriqs.gf
import pytriqs.lattice

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/types.hpp") # Manually added
module.add_include("../c++/lattice.hpp") # Manually added

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
#include <triqs/cpp2py_converters/variant.hpp>

using namespace triqs::gfs;
using namespace triqs::lattice;
using namespace tprf;
""")

module.add_function("gk_iw_t g0k_from_ek(double mu, gf_view<brillouin_zone, matrix_valued> ek, gf_mesh<imfreq> mesh)", doc = """""")
module.add_function("gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, gf_view<imfreq, matrix_valued> sigma)", doc = """""")

module.add_function("gr_iw_t gr_from_gk(gf_view<cartesian_product<imfreq, brillouin_zone>> gk)", doc = """""")
module.add_function("gk_iw_t gk_from_gr(gf_view<cartesian_product<imfreq, cyclic_lattice>> gr, brillouin_zone bz)", doc = """""")

module.add_function("chi0r_t chi0r_from_gr_PH(int nw, int nnu, gf_view<cartesian_product<imfreq, cyclic_lattice>> gr)", doc = """""")
module.add_function("chi0r_t chi0r_from_chi0q(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")
module.add_function("chi0q_t chi0q_from_chi0r(gf_view<cartesian_product<imfreq, imfreq, cyclic_lattice>, tensor_valued<4>> chi0r, brillouin_zone bz)", doc = """""")

module.add_function ("gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")

module.generate_code()
