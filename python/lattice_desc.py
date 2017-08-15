# Generated automatically using the command :
# c++2py.py -mpytriqs.applications.tprf.lattice ../c++/lattice.hpp -I../c++ --moduledoc "Green's function lattice tools"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.tprf.lattice", doc = "Green's function lattice tools", app_name = "pytriqs.applications.tprf.lattice")

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


# cpp2py does not cope with typedefs... at all
#module.add_function ("gk_iw_t g0k_from_ek(double mu, ek_vt ek, g_iw_t::mesh_t mesh)", doc = """""")
#module.add_function ("gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, g_iw_vt sigma)", doc = """""") 

module.add_function("gk_iw_t g0k_from_ek(double mu, gf_view<brillouin_zone, matrix_valued> ek, gf_mesh<imfreq> mesh)", doc = """""")
module.add_function("gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, gf_view<imfreq, matrix_valued> sigma)", doc = """""")

module.add_function("gr_iw_t gr_from_gk(gf_view<cartesian_product<imfreq, brillouin_zone>> gk)", doc = """""")
module.add_function("gk_iw_t gk_from_gr(gf_view<cartesian_product<imfreq, cyclic_lattice>> gr, brillouin_zone bz)", doc = """""")

module.add_function("chi0r_t chi0r_from_gr_PH(int nw, int nnu, gf_view<cartesian_product<imfreq, cyclic_lattice>> gr)", doc = """""")
module.add_function("chi0r_t chi0r_from_chi0q(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")
module.add_function("chi0q_t chi0q_from_chi0r(gf_view<cartesian_product<imfreq, imfreq, cyclic_lattice>, tensor_valued<4>> chi0r, brillouin_zone bz)", doc = """""")

module.add_function ("chi0_t get_at_q(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q, mini_vector<int, 3> q)", doc = """""")
module.add_function ("gf<cartesian_product<imfreq>, tensor_valued<4>> chi0_sum_nu(gf_view<cartesian_product<imfreq, imfreq>, tensor_valued<4>> chi0)", doc = """""")

module.add_function ("gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")

#module.add_function ("", doc = """""")
#module.add_function ("", doc = """""")
#module.add_function ("", doc = """""")
#module.add_function ("", doc = """""")
#module.add_function ("", doc = """""")
#module.add_function ("", doc = """""")
#module.add_function ("", doc = """""")

module.generate_code()
