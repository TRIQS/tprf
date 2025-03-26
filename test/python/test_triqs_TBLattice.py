################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2025 by Xaver Landerl
# Author: Xaver Landerl
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

"""
Test Triqs tight binding object comparing
dispersions for square and triangular lattice
with respect to analytical expressions.
"""

import numpy as np

from triqs.lattice.tight_binding import TBLattice


def e_k_square(kx, ky):
    return -2*(np.cos(kx) + np.cos(ky))


def e_k_triangular(kx, ky):
    return -2*(np.cos(kx) + 2*np.cos(kx/2) * np.cos(np.sqrt(3)/2*ky))


def tight_binding_square():
    units = [(1,0,0), (0,1,0)]
    hoppings = {(+1,+0) : [[-1.0]],
                (-1,+0) : [[-1.0]],
                (+0,+1) : [[-1.0]],
                (+0,-1) : [[-1.0]]}
    Lat = TBLattice(units=units, hoppings=hoppings)
    return Lat


def tight_binding_triangular():
    units = [(1,0,0), (1/2,np.sqrt(3)/2,0)]
    hoppings = {(+1,+0) : [[-1.0]],
                (-1,+0) : [[-1.0]],
                (+0,+1) : [[-1.0]],
                (+0,-1) : [[-1.0]],
                (+1,-1) : [[-1.0]],
                (-1,+1) : [[-1.0]]}
    Lat = TBLattice(units=units, hoppings=hoppings)
    return Lat


def run_dispersion(disp='square', verbose=False):

    if disp == 'square':
        e_k_anal = e_k_square
        Lat = tight_binding_square()
        
    elif disp == 'triangular':
        e_k_anal = e_k_triangular
        Lat = tight_binding_triangular()

    else:
        raise NotImplementedError
    
    k_mesh = Lat.get_kmesh(n_k=512)
    e_k = Lat.fourier(k_mesh)

    e_k_interpolator = np.vectorize(lambda kx, ky : float(e_k((kx, ky, 0))[0, 0].real))
    
    k = np.linspace(-np.pi, np.pi, 100)
    KX, KY = np.meshgrid(k, k)

    e_k_interp = e_k_interpolator(KX, KY)
    e_k_ref = e_k_anal(KX, KY)

    if verbose:
        print(np.max(np.abs(e_k_interp - e_k_ref)))
        print(np.linalg.norm(e_k_interp - e_k_ref))

    np.testing.assert_array_almost_equal(e_k_interp, e_k_ref, decimal=3)


def test_square_dispersion(verbose=False):
    run_dispersion(disp='square')
    

def test_triangular_dispersion(verbose=False):
    run_dispersion(disp='triangular')

    
if __name__ == '__main__':
    test_square_dispersion()
    test_triangular_dispersion()
