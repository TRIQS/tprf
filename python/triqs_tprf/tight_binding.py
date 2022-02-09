################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
# Copyright (C) 2018 The Simons Foundation
# Copyright (C) 2019, S. Käser
# Author: Hugo U. R. Strand, S. Käser
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import itertools

from triqs.lattice.tight_binding import TBLattice


# ----------------------------------------------------------------------
def create_square_lattice(norb, t, tp=0.0, zeeman=0.0, spin=False, **kwargs):
    r"""Retuns TBLattice that represents a model on a square lattice
    
    The model is described by the Hamiltonian

    .. math:: 
        H=-t \sum_{\langle j, l\rangle \sigma}\left(c_{j \sigma}^{\dagger}
        c_{l \sigma}+c_{l \sigma}^{\dagger} c_{j \sigma}\right) -
        t' \sum_{\langle\langle j, l\rangle\rangle \sigma}\left(c_{j \sigma}^{\dagger}
        c_{l \sigma}+c_{l \sigma}^{\dagger} c_{j \sigma}\right) +
        \xi \sum_j \left(n_{j \uparrow} - n_{j \downarrow} \right)\,,

    where the angular bracket describes a sum over nearest neighbors and the double
    angular bracket over next-nearest neighbors.


    Parameters
    ----------
    norb : int,
           Number of orbitals excluding spin
    t : complex,
        Kinetic energy of nearest neighbor hopping.
        Corresponds to :math:`t`.
    tp : complex, optional
         Kinetic energy of next-nearest neighbor hopping.
         Corresponds to :math:`t'`.
    zeeman : complex, optional
             Strength of Zeeman term.
             Corresponds to :math:`\xi`.
    spin : bool,
           True if spin index should be used explicitly, False otherwise.
           The Zeeman term can only be applied if spin=True.

    Returns
    -------
    square_lattice : TBLattice
    """

    if zeeman != 0.0 and not spin:
        raise AttributeError('There can not be a Zeeman term in a spinless model.')
    if spin:
        norb *= 2

    import numpy as np
    t_matrix = -t * np.eye(norb)
    tp_matrix = -tp * np.eye(norb)
    zeeman_matrix = zeeman * np.diag([(-1)**orb for orb in range(norb)])

    hopping = {
                # Zeeman term
                ( 0, 0): zeeman_matrix,

                # nearest neighbour hopping
                ( 0,+1): t_matrix,
                ( 0,-1): t_matrix,
                (+1, 0): t_matrix,
                (-1, 0): t_matrix,
                
                # next-nearest neighbour hopping
                ( +1,+1): tp_matrix,
                ( -1,-1): tp_matrix,
                (+1, -1): tp_matrix,
                (-1, +1): tp_matrix,
                }

    units = [(1, 0, 0), (0, 1, 0)]
    orbital_positions = [(0, 0, 0)] * norb

    square_lattice = TBLattice(units, hopping, orbital_positions)

    return square_lattice

# ----------------------------------------------------------------------
def create_model_for_tests(norb, dim, t=0, t1=0, t2=0, t12=0, t21=0, **kwargs):
    full_units = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    units = full_units[:dim]

    import numpy as np
    if norb == 1:
        t_matrix = -t * np.eye(norb)
    elif norb == 2:
        t_matrix = -np.array([[t1, t12], [t21, t2]])
    else:
        raise NotImplementedError

    all_nn_hoppings = list(itertools.product([-1, 0, 1], repeat=dim)) 
    non_diagonal_hoppings = [hopping for hopping in all_nn_hoppings if sum(np.abs(hopping)) == 1] 
    hoppings = {hopping : t_matrix for hopping in non_diagonal_hoppings}

    orbital_positions = [(0, 0, 0)] * norb

    model_for_tests = TBLattice(units, hoppings, orbital_positions)

    return model_for_tests

__all__ = ['TBLattice', 'create_square_lattice', 'create_model_for_tests']
