# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/hubbard_atom.hpp --members_read_only -N triqs_tprf::hubbard_atom -a triqs_tprf -m hubbard_atom -o hubbard_atom -C pytriqs --moduledoc="Exact correlation functions for the hubbard atom" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "hubbard_atom", doc = "Exact correlation functions for the hubbard atom", app_name = "triqs_tprf")

# Imports
import pytriqs.gf

# Add here all includes
module.add_include("triqs_tprf/hubbard_atom.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf::hubbard_atom;
""")


module.add_function ("triqs_tprf::g_iw_t triqs_tprf::hubbard_atom::single_particle_greens_function (int nw, double beta, double U)", doc = """Single-particle Green\'s function of the Hubbard atom at half-filling\n\n     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)\n     please cite the paper if you use this function!\n\n     .. math::\n         G(i\\omega_n) = \\frac{1}{i\\omega_n - \\frac{U^2}{4 i\\omega_n}}\n\n     :param nw: number of Matsubara frequencies\n     :param beta: inverse temperature\n     :param U: Hubbard U interaction parmeter\n     @return g single-particle Green\'s function of the Hubbard atom :math:`G(i\\omega_n)`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::hubbard_atom::chi_ph_magnetic (int nw, int nwf, double beta, double U)", doc = """Magnetic susceptibility of the Hubbard atom at half-filling :math:`\\chi(\\omega, \\nu, \\nu\')`\n\n     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)\n     please cite the paper if you use this function!\n\n     :param nw: number of bosonic Matsubara frequencies\n     :param nwf: number of fermionic Matsubara frequencies\n     :param beta: inverse temperature\n     :param U: Hubbard U interaction parmeter\n     @return chi magnetic susceptibility :math:`\\chi(\\omega, \\nu, \\nu\')`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::hubbard_atom::gamma_ph_magnetic (int nw, int nwf, double beta, double U)", doc = """Magnetic vertex function in the particle-hole channel of the Hubbard atom at half-filling :math:`\\Gamma(\\omega, \\nu, \\nu\')`\n\n    Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)\n    please cite the paper if you use this function!\n\n    :param nw: number of bosonic Matsubara frequencies\n    :param nwf: number of fermionic Matsubara frequencies\n    :param beta: inverse temperature\n    :param U: Hubbard U interaction parmeter\n    @return gamma magnetic susceptibility :math:`\\Gamma(\\omega, \\nu, \\nu\')`""")



module.generate_code()