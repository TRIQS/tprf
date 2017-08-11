# Generated automatically using the command :
# c++2py.py -mpytriqs.applications.tprf.fourier ../c++/fourier.hpp -I../c++ --moduledoc "Green's function FFT tools"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.tprf.fourier", doc = "Green's function FFT tools", app_name = "pytriqs.applications.tprf.fourier")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/fourier.hpp") # Manually added

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/gf.hpp>
#include <triqs/python_tools/converters/block_gf.hpp>
#include <triqs/python_tools/converters/block2_gf.hpp>
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/set.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/optional.hpp>
#include <triqs/python_tools/converters/variant.hpp>
#include <triqs/python_tools/converters/tuple.hpp>
using namespace tprf::fourier;
using namespace tprf;
""")

module.add_function ("chi4_iw_t chi4_iw_from_tau (gf_view<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>> chi4_tau, int n_iw)", doc = """""")
module.add_function ("chi3_iw_t chi3_iw_from_tau (gf_view<cartesian_product<imtime, imtime>, tensor_valued<4>> chi3_tau, int n_iw)", doc = """""")
module.add_function ("chi2_iw_t chi2_iw_from_tau (gf_view<imtime, tensor_valued<4>> chi2_tau, int n_iw)", doc = """""")
module.add_function ("g_iw_t g_iw_from_tau (gf_view<imtime, matrix_valued> g_tau, int n_iw)", doc = """""")

module.add_function ("chi4_tau_t chi4_tau_from_iw (gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> chi4_iw, int n_tau)", doc = """""")
module.add_function ("chi3_tau_t chi3_tau_from_iw (gf_view<cartesian_product<imfreq, imfreq>, tensor_valued<4>> chi3_iw, int n_tau)", doc = """""")
module.add_function ("chi2_tau_t chi2_tau_from_iw (gf_view<imfreq, tensor_valued<4>> chi2_iw, int n_tau)", doc = """""")
module.add_function ("g_tau_t g_tau_from_iw (gf_view<imfreq, matrix_valued> g_iw, int n_tau)", doc = """""")

module.generate_code()
