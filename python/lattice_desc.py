
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

module.add_function("array<std::complex<double>, 6> cluster_mesh_fourier_interpolation(array<double, 2> k_vecs, chi_wr_cvt chi)", doc = """""")


module.add_function("gk_iw_t lattice_dyson_g0_wk(double mu, ek_vt e_k, g_iw_t::mesh_t mesh)", doc = """""")
module.add_function("gk_iw_t lattice_dyson_g_wk(double mu, ek_vt e_k, gf_view<imfreq, matrix_valued> sigma_w)", doc = """""")
module.add_function("g_iw_t lattice_dyson_g_w(double mu, ek_vt e_k, g_iw_vt sigma_w)", doc = """""")


module.add_function("gr_iw_t fourier_wk_to_wr(gf_view<cartesian_product<imfreq, brillouin_zone>> g_k)", doc = """""")
module.add_function("gk_iw_t fourier_wr_to_wk(gf_view<cartesian_product<imfreq, cyclic_lattice>> g_r)", doc = """""")
module.add_function("gr_tau_t fourier_wr_to_tr(gr_iw_vt g_wr, int ntau=-1)", doc = """""")


# -- Bubble in imaginary time

module.add_function("chi_tr_t chi0_tr_from_grt_PH(gr_tau_vt grt)", doc = """""")
module.add_function("chi_wr_t chi0_w0r_from_grt_PH(gr_tau_vt grt)", doc = """""")
module.add_function("chi_wr_t chi_w0r_from_chi_tr(chi_tr_vt chi_tr)", doc = """""")
module.add_function("chi_wr_t chi_wr_from_chi_tr(chi_tr_vt chi_tr, int nw)", doc = """""")
module.add_function("chi_wk_t chi_wk_from_chi_wr(chi_wr_vt chi_wr)", doc = """""")
module.add_function("chi_wr_t chi_wr_from_chi_wk(chi_wk_vt chi_wk)", doc = """""")

# -- Bubble static analytic

module.add_function("chi_wk_t lindhard_chi00_wk(gf<brillouin_zone, matrix_valued> e_k, int nw, double beta, double mu)", doc = """""")

# -- RPA

module.add_function("chi_wk_t solve_rpa_PH(chi_wk_vt chi0, array_view<std::complex<double>, 4> U)", doc = """""")

# -- Full BSE functions

module.add_function("chi0r_t chi0r_from_gr_PH(int nw, int nnu, gf_view<cartesian_product<imfreq, cyclic_lattice>> gr)", doc = """""")
module.add_function("chi0q_t chi0q_from_g_wk_PH(int nw, int nnu, gf_view<cartesian_product<imfreq, brillouin_zone>> g_wk)", doc = """""")

module.add_function("chi0r_t chi0r_from_chi0q(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")
module.add_function("chi0q_t chi0q_from_chi0r(gf_view<cartesian_product<imfreq, imfreq, cyclic_lattice>, tensor_valued<4>> chi0r)", doc = """""")

module.add_function ("gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")
module.add_function ("gf<cartesian_product<imfreq, brillouin_zone>, tensor_valued<4>> chi0q_sum_nu_tail_corr_PH(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q)", doc = """""")

module.add_function ("gf<imfreq, tensor_valued<4>> chi0q_sum_nu_q(chi0q_t chi0q)", doc = """""")

module.add_function("chiq_t chiq_from_chi0q_and_gamma_PH(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> gamma_ph)", doc = """""")

module.add_function("gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu_from_chi0q_and_gamma_PH(gf_view<cartesian_product<imfreq, imfreq, brillouin_zone>, tensor_valued<4>> chi0q, gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> gamma_ph)", doc = """""")

module.add_function("gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu_from_g_wk_and_gamma_PH(gk_iw_t g_wk, g2_iw_vt gamma_ph_wnn, int tail_corr_nwf=-1)", doc = """""")

module.add_function("gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH(double mu, ek_vt e_k, g_iw_vt sigma_w, g2_iw_vt gamma_ph_wnn, int tail_corr_nwf=-1)", doc = """""")

module.add_function ("gf<cartesian_product<brillouin_zone, imfreq>, tensor_valued<4>> chiq_sum_nu(gf_view<cartesian_product<brillouin_zone, imfreq, imfreq, imfreq>, tensor_valued<4>> chiq)", doc = """""")

module.add_function ("gf<imfreq, tensor_valued<4>> chiq_sum_nu_q(gf_view<cartesian_product<brillouin_zone, imfreq, imfreq, imfreq>, tensor_valued<4>> chiq)", doc = """""")

module.generate_code()
