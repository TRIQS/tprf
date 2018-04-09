
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
#include <cpp2py/converters/variant.hpp>

using namespace triqs::gfs;
using namespace triqs::lattice;
using namespace tprf;
""")

module.add_function("gk_iw_t g0k_from_ek(double mu, gf_view<brillouin_zone, matrix_valued> ek, gf_mesh<imfreq> mesh)", doc = """""")
module.add_function("gk_iw_t gk_from_ek_sigma(double mu, ek_vt ek, gf_view<imfreq, matrix_valued> sigma)", doc = """""")

module.add_function("gr_iw_t gr_from_gk(gf_view<cartesian_product<imfreq, brillouin_zone>> gk)", doc = """""")
module.add_function("gk_iw_t gk_from_gr(gf_view<cartesian_product<imfreq, cyclic_lattice>> gr)", doc = """""")

# -- Bubble in imaginary time

module.add_function("gr_tau_t grt_from_grw(gf_view<cartesian_product<imfreq, cyclic_lattice>> grw)", doc = """""")
module.add_function("chi_tr_t chi0_tr_from_grt_PH(gr_tau_vt grt)", doc = """""")
module.add_function("chi_wr_t chi_wr_from_chi_tr(chi_tr_vt chi_tr)", doc = """""")
module.add_function("chi_wk_t chi_wk_from_chi_wr(chi_wr_vt chi_wr)", doc = """""")

# -- Bubble static analytic

module.add_function("chi_wk_t chi00_wk_from_ek(gf<brillouin_zone, matrix_valued> ek_in, int nw, double beta, double mu)", doc = """""")

# -- Full BSE functions

module.add_function("chi0r_t chi0r_from_gr_PH(int nw, int nnu, gf_view<cartesian_product<imfreq, cyclic_lattice>> gr)", doc = """""")
module.add_function("chi0r_t chi0r_from_chi0q(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")
module.add_function("chi0q_t chi0q_from_chi0r(gf_view<cartesian_product<imfreq, imfreq, cyclic_lattice>, tensor_valued<4>> chi0r)", doc = """""")

module.add_function ("gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")
module.add_function ("gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu_tail_corr_PH(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")

module.add_function ("gf<imfreq, tensor_valued<4>> chi0q_sum_nu_q(chi0q_t chi0q)", doc = """""")

module.add_function("chiq_t chiq_from_chi0q_and_gamma_PH(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> gamma_ph)", doc = """""")

module.add_function ("gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu(gf_view<cartesian_product<brillouin_zone, imfreq, imfreq, imfreq>, tensor_valued<4>> chiq)", doc = """""")

module.add_function ("gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(gf_view<cartesian_product<brillouin_zone, imfreq, imfreq, imfreq>, tensor_valued<4>> chiq)", doc = """""")

module.generate_code()
