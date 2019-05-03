# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/hubbard_atom.hpp --members_read_only -N triqs_tprf::hubbard_atom -a triqs_tprf -m hubbard_atom -o hubbard_atom -C pytriqs --moduledoc="Exact correlation functions for the hubbard atom" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "hubbard_atom", doc = r"Exact correlation functions for the hubbard atom", app_name = "triqs_tprf")

# Imports
module.add_imports(*['pytriqs.gf'])

# Add here all includes
module.add_include("triqs_tprf/hubbard_atom.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf::hubbard_atom;
""")


module.add_function ("triqs_tprf::g_iw_t triqs_tprf::hubbard_atom::single_particle_greens_function (int nw, double beta, double U)", doc = r"""Single-particle Green's function of the Hubbard atom at half-filling

     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
     please cite the paper if you use this function!

     .. math::
         G(i\omega_n) = \frac{1}{i\omega_n - \frac{U^2}{4 i\omega_n}}

Parameters
----------
nw
     number of Matsubara frequencies

beta
     inverse temperature

U
     Hubbard U interaction parmeter

Returns
-------
out
     g single-particle Green's function of the Hubbard atom :math:`G(i\omega_n)`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::hubbard_atom::chi_ph_magnetic (int nw, int nwf, double beta, double U)", doc = r"""Magnetic susceptibility of the Hubbard atom at half-filling :math:`\chi(\omega, \nu, \nu')`

     Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
     please cite the paper if you use this function!

Parameters
----------
nw
     number of bosonic Matsubara frequencies

nwf
     number of fermionic Matsubara frequencies

beta
     inverse temperature

U
     Hubbard U interaction parmeter

Returns
-------
out
     chi magnetic susceptibility :math:`\chi(\omega, \nu, \nu')`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::hubbard_atom::gamma_ph_magnetic (int nw, int nwf, double beta, double U)", doc = r"""Magnetic vertex function in the particle-hole channel of the Hubbard atom at half-filling :math:`\Gamma(\omega, \nu, \nu')`

    Using analytical formulas from Thunstrom et al. PRB 98, 235107 (2018)
    please cite the paper if you use this function!

Parameters
----------
nw
     number of bosonic Matsubara frequencies

nwf
     number of fermionic Matsubara frequencies

beta
     inverse temperature

U
     Hubbard U interaction parmeter

Returns
-------
out
     gamma magnetic susceptibility :math:`\Gamma(\omega, \nu, \nu')`""")



module.generate_code()