################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
# Copyright (C) 2018 The Simons Foundation
# Author: Hugo U. R. Strand
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

import numpy as np

from pytriqs.lattice.lattice_tools import BrillouinZone as BrillouinZone
from pytriqs.lattice.lattice_tools import BravaisLattice as BravaisLattice
from pytriqs.lattice.lattice_tools import TightBinding as TightBinding

from pytriqs.lattice.lattice_tools import dos as dos_c
from pytriqs.lattice.lattice_tools import dos_patch as dos_patch_c

from pytriqs.lattice.lattice_tools import energies_on_bz_grid, energies_on_bz_path, hopping_stack, energy_matrix_on_bz_path

from pytriqs.dos import DOS
from pytriqs.gf import Gf, MeshBrillouinZone, MeshCyclicLattice

class TBLattice(object):

    """ Class describing a tight binding lattice model.

    This class is based in the TRIQS tight binding class and has been extended with
    some extra convienience methods and documentation.

    Parameters
    ----------

    units : list of three-tuples of floats
        Basis vectors for the real space lattice.

    hopping : dict
        Dictionary with three tuple of integeers as keys, 
        describing real space hoppings in multiples of 
        the real space lattice basis vectors, and values being 
        numpy ndarray hopping matrices in the orbital indices.

    orbital_positions : list of three three-tuples of floats.
        Internal orbital positions in the unit-cell.

    orbital_names : list of strings
        Names for each orbital.

    """

    def __init__ (self, units, hopping, orbital_positions = [ (0, 0, 0) ], orbital_names = [""]):

        # the k are int32 which boost python does like to convert 
        def reg(k) : return tuple( int(x) for x in k) 
        self._hop = dict ( ( reg(k), np.array(v)) for k, v in hopping.items())
        orb = dict ( (str(i), orb) for (i, orb) in enumerate(orbital_positions ))
        self.bl = BravaisLattice(units, orbital_positions)
        self.bz = BrillouinZone(self.bl)
        self.tb = TightBinding(self.bl, self._hop) #, orbital_positions )
        self.dim = self.bl.dim
        self.NOrbitalsInUnitCell = self.bl.n_orbitals
        self.Units = units
        self.OrbitalPositions = orbital_positions 
        self.OrbitalNames = orbital_names
        self.MuPattern = np.identity(self.NOrbitalsInUnitCell)

    def latt_to_real_x(self, p) : 
        return self.bl.lattice_to_real_coordinates (np.array(p, np.float64))

    def hopping_dict(self) : return self._hop

    def hopping(self, k_stack) :
        return hopping_stack(self.tb, k_stack)

    def periodization_matrix(self, n_k):
        n_k = np.array(n_k)
        assert( len(n_k) == 3 )
        assert( n_k.dtype == np.int )
        periodization_matrix = np.diag(np.array(list(n_k), dtype=np.int32))
        return periodization_matrix
    
    def get_kmesh(self, n_k):
        return MeshBrillouinZone(self.bz, self.periodization_matrix(n_k))

    def get_rmesh(self, n_k):
        return MeshCyclicLattice(self.bl, self.periodization_matrix(n_k))
    
    def on_mesh_brillouin_zone(self, n_k):

        """ Construct a discretization of the tight binding model in
        reciprocal space with given number of k-points.

        Parameters
        ----------

        n_k : three-tuple of ints
            Number of k-points in every dimension.

        Returns
        -------

        e_k : TRIQS Greens function on a Brillioun zone mesh
            Reciprocal space tight binding dispersion.

        """
        
        target_shape = [self.NOrbitalsInUnitCell] * 2

        kmesh = self.get_kmesh(n_k)

        e_k = Gf(mesh=kmesh, target_shape=target_shape)

        k_vec = np.array([k.value for k in kmesh])
        k_vec_rel = np.dot(np.linalg.inv(self.bz.units()).T, k_vec.T).T   
        e_k.data[:] = self.hopping(k_vec_rel.T).transpose(2, 0, 1)

        return e_k

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
