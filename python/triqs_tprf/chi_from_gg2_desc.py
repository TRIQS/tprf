# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/chi_from_gg2.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m chi_from_gg2 -o chi_from_gg2 -C pytriqs --moduledoc="Calculation of generalized susceptibilities" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "chi_from_gg2", doc = r"Calculation of generalized susceptibilities", app_name = "triqs_tprf")

# Imports
module.add_imports(*['pytriqs.gf'])

# Add here all includes
module.add_include("triqs_tprf/chi_from_gg2.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", doc = r"""Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::chi0_from_gg2_PH (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2)", doc = r"""Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Hole channel

     Computes

     .. math::
         \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega + \nu)

Parameters
----------
g
     single particle Green's function :math:`G_{ab}(\nu)`

g2
     two-particle Green's function with the mesh to use for
     :math:`\chi^{(0)}`

Returns
-------
out
     chi0 particle-hole bubble
     :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::chi0_from_gg2_PP (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2)", doc = r"""Bubble susceptibility :math:`\chi^{(0)} = GG` in the Particle-Particle
   channel

     Computes

     .. math::
         \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         - \beta \delta_{\nu, \nu'} G_{da}(\nu) \cdot G_{bc}(\omega - \nu)

Parameters
----------
g
     single particle Green's function :math:`G_{ab}(\nu)`

g2
     two-particle Green's function with the mesh to use for
     :math:`\chi^{(0)}`

Returns
-------
out
     chi0 particle-particle bubble
     :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::chi_from_gg2_PH (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2)", doc = r"""Generalized susceptibility :math:`\chi^{(0)} = G^{(2)} - GG` in the
   Particle-Hole channel

     Computes

     .. math::
         \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
         - \beta \delta_{\omega} G_{ba}(\nu) \cdot G_{dc}(\nu')

Parameters
----------
g
     single particle Green's function :math:`G_{ab}(\nu)`

g2
     two-particle Green's function
     :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`

Returns
-------
out
     chi generalized particle-hole susceptibility
     :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`""")

module.add_function ("triqs_tprf::g2_iw_t triqs_tprf::chi_from_gg2_PP (triqs_tprf::g_iw_vt g, triqs_tprf::g2_iw_vt g2)", doc = r"""Generalized susceptibility :math:`\chi^{(0)} = G^{(2)} - GG` in the
   Particle-Particle channel

     Computes

     .. math::
         \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu') =
         G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')
         - \beta \delta_{\nu + \nu' - \omega} G_{ba}(\nu) \cdot G_{dc}(\nu')

Parameters
----------
g
     single particle Green's function :math:`G_{ab}(\nu)`

g2
     two-particle Green's function
     :math:`G^{(2)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`

Returns
-------
out
     chi generalized particle-hole susceptibility
     :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu,\nu')`""")



module.generate_code()