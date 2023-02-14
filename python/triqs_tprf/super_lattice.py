
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

import numpy

#from triqs.lattice.tight_binding import TBLattice
from triqs_tprf.tight_binding import TBLattice

__all__ = ['TBSuperLattice']


def nint_strict(x, precision=1.e-9):
    """ Round array x to the closest integer and asserts that its distance to this integer is less than precision.
        precision must satisfy :  precision >0 and precision <0.5
    """
    assert precision >0 and precision <0.5, "nint_strict : precision makes no sense !"
    i = numpy.floor(x+0.5)
    assert abs(i-x).max() < precision, repr(i) + "\n "+repr(x) + "\n The Float is not close enough to the integer "
    return i.astype(int)


class TBSuperLattice(TBLattice):
    
    r""" Class building a superlattice on top of a base lattice (TBLattice).

    The class inherits from TBLattice and has all the basic TBLattice methods, especially the ``on_mesh_brillouin_zone``-method.

    Parameters
    ----------

    tb_lattice : TBLattice instance
        The base tight binding lattice.

    super_lattice_units : ndarray (2D)
        The unit vectors of the superlattice in the ``tb_lattice`` (integer) coordinates.

    cluster_sites : 
        Coordinates of the cluster in tb_lattice coordinates. 
        If ``None``, an automatic computation of cluster positions 
        is made as follows: it takes all points whose coordinates 
        in the basis of the superlattice are in [0, 1[^dimension.

    remove_internal_hoppings : bool
        If ``true``, the hopping terms are removed inside the cluster. 
        Useful to add Hartree Fock terms at the boundary of a cluster, e.g.

    """    
    def __init__(self, tb_lattice, super_lattice_units, cluster_sites = None, remove_internal_hoppings = False):

        #if not isinstance(tb_lattice, TBLattice): raise ValueError, "tb_lattice should be an instance of TBLattice"
        self.__BaseLattice = tb_lattice
        dim = tb_lattice.dim

        try:
            self.__super_lattice_units = numpy.array(super_lattice_units, copy=True)
            assert self.__super_lattice_units.shape == (dim, dim)
        except:
            raise ValueError("super_lattice_units is not correct. Cf Doc. value is %s, dim = %s "%(super_lattice_units,dim))

        Ncluster_sites = int(numpy.rint(abs(numpy.linalg.det(self.__super_lattice_units ))))
        assert Ncluster_sites >0, "Superlattice vectors are not independant !"
        self._M = self.__super_lattice_units.transpose()
        self._Mtilde = numpy.array(numpy.rint(numpy.linalg.inv(self._M)*Ncluster_sites), dtype = int)

        self.__remove_internal_hoppings = remove_internal_hoppings
        #self.norb = tb_lattice.NOrbitalsInUnitCell
        self.Norb = tb_lattice.NOrbitalsInUnitCell * Ncluster_sites

        # cluster_sites computation
        if cluster_sites!=None:
            self.__cluster_sites = list(cluster_sites)[:]
        else: # Computes the position of the cluster automatically
            self.__cluster_sites = []
            #We tile the super-cell with the tb_lattice points and retains
            # the points inside it and store it.
            #M=numpy.array(self.__super_lattice_units) # BUG!
            M=numpy.array(self.__super_lattice_units.transpose())
            assert M.shape==tuple(2*[dim]), "Tiling Construct: super_lattice_units does not have the correct size"
            #Minv = Ncluster_sites*numpy.linalg.inverse(M)  #+0.5  # +0.5 is for the round
            #Mtilde = Minv.astype(numpy.Int)  # now is integer.
            Mtilde = nint_strict(Ncluster_sites*numpy.linalg.inv(M))
            #print 'Mtilde (inside cluster sites) =\n', Mtilde.__repr__()
            # round to the closest integer, with assert that precision is <1.e-9
            if dim==1:  a=(max(M[0,:]), 0, 0 )
            elif dim==2:  a=(2*max(M[0,:]), 2*max(M[1,:]), 0 )
            elif dim==3: a= (3*max(M[0,:]), 3*max(M[1,:]), 3*max(M[2,:]))
            else: raise ValueError("dim is not between 1 and 3 !!")
            r = lambda i:  list(range(-a[i] , a[i]+1))
            for nx in r(0):
                for ny in r(1):
                    for nz in r(2):
                        nv = numpy.array([nx, ny, nz][0:dim])
                        k_sl = numpy.dot(Mtilde, nv)
                        if (min(k_sl) >= 0) and (max(k_sl) < Ncluster_sites ): # The point is a point of the cluster. We store it.
                            self.__cluster_sites.append(nv.tolist())

        assert len(self.__cluster_sites) == Ncluster_sites, """Number of cluster positions incorrect (compared to the volume of unit cell of the Superlattice)"""
        self.Ncluster_sites = Ncluster_sites

        # creating a dictionnary position_of_sites -> number e.g. (1, 0): 2 etc...
        # self._clustersites_hash =  dict ([ (tuple(pos), n) for n, pos in enumerate(self.cluster_sites)])

        #print 'Ns = ', self.Ncluster_sites
        #print 'cluster_sites =', self.__cluster_sites
        #print 'M =\n', self._M.__repr__()
        #print 'Mtilde =\n', self._Mtilde.__repr__()
        #import numpy as np
        #print 'M*Mtilde =\n', np.dot(self._M, self._Mtilde)
        #exit()
        
        # Compute the new Hopping in the supercell
        Hopping = self.fold(tb_lattice.hopping_dict(), remove_internal_hoppings)
        if 0:
            for k, v in list(Hopping.items()):
                print(k)
                print(v.real)

        # Compute the new units of the lattice in real coordinates
        Units = numpy.dot(self.__super_lattice_units, tb_lattice.Units)

        # Positions and names of orbitals in the supercell: just translate all orbitals for cluster site positions
        # in R^3 coordinates.
        Orbital_Positions = [POS + tb_lattice.latt_to_real_x(CS) for POS in tb_lattice.OrbitalPositions for CS in self.__cluster_sites]

        #Orbital_Names = [ '%s%s'%(n, s) for n in tb_lattice.OrbitalNames for s in range(Ncluster_sites)]
        site_index_list, orbital_index_list = list(range(1, Ncluster_sites+1)), tb_lattice.OrbitalNames
        if len(orbital_index_list)==1:
            Orbital_Names= [ s for s in site_index_list ]
        elif len(site_index_list)==1 and len(orbital_index_list)>1:
            Orbital_Names= [ o for o in orbital_index_list]
        elif len(site_index_list)>1 and len(orbital_index_list)>1:
            Orbital_Names= [ (pos, o) for o in orbital_index_list for pos in site_index_list]

        #print tb_lattice.OrbitalNames #Orbital_Names
        TBLattice.__init__(self, Units, Hopping, Orbital_Positions, Orbital_Names)
        # we pass False since the folding has arealdy been done in tb_lattice

        assert self.Norb == self.NOrbitalsInUnitCell

    __HDF_reduction__ = ['__BaseLattice', '__super_lattice_units', '__cluster_sites', '__remove_internal_hoppings']

    def __reduce__ (self):
        return tuple([getattr(self, x) for x in self.__HDF_reduction__])

    def fold(self, D1, remove_internal=False, create_zero = None):
        """ Input: a function  r-> f(r) on the tb_lattice given as a dictionnary
            Output: the function R-> F(R) folded on the superlattice.
            Only requirement is that f(r)[orbital1, orbital2] is properly defined.
            Hence f(r) can be a numpy, a GFBloc, etc...
            """
        #Res , norb = {} , self.__BaseLattice.NOrbitalsInUnitCell
        Res , norb = {} , len(list(D1.values())[0])
        pack = self.pack_index_site_orbital
        for nsite, CS in enumerate(self.__cluster_sites):
            for disp, t in list(D1.items()):
                #print 'CS, disp =', CS, disp
                R, alpha = self.change_coordinates_L_to_SL(numpy.array(CS)+numpy.array(disp))
                if R not in Res: Res[R] = create_zero() if create_zero else numpy.zeros((self.Norb, self.Norb), dtype = type(t[0,0]))
                if not(remove_internal) or R!= self.tb_lattice.dim*(0, ):
                    for orb1 in range(norb):
                        for orb2 in range(norb):
                            Res[R][pack(nsite, orb1), pack(alpha, orb2)] += t[orb1, orb2]
        return Res

    def change_coordinates_SL_to_L(self, R , alpha):
        """Given a point in the supercell R, site (number) alpha, it computes its position on the tb_lattice in lattice coordinates"""
        return numpy.dot (self._M, numpy.array(R)) + self.__cluster_sites[alpha,:]

    def change_coordinates_L_to_SL(self, x):
        """Given a point on the tb_lattice in lattice coordinates, returns its coordinates (R, alpha) in the Superlattice"""
        aux  = numpy.dot(self._Mtilde, numpy.array(x))
        R = aux // self.Ncluster_sites
        dx = list (x - numpy.dot (self._M, R) ) # force int ?

        if False:
            #dr = numpy.dot(self._Mtilde, dx)
            #dr = aux - self.Ncluster_sites * R

            print('M * R =', numpy.dot(self._M, R))
            print('aux =', aux)
            print('R, dx =', R, dx)
            print('dr =', dr)
        
        return tuple(R), self.__cluster_sites.index(dx)

    def pack_index_site_orbital(self, n_site, n_orbital):
        """ nsite and n_orbital must start at 0"""
        return n_site + (n_orbital ) * self.Ncluster_sites

    def unpack_index_site_orbital (self, index):
        """Inverse of pack_index_site_orbital"""
        n_orbital  =   (index)//self.Ncluster_sites
        n_site =  index - n_orbital*self.Ncluster_sites
        return n_site, n_orbital

    def cluster_sites(self):
        """
           Generate the position of the cluster site in the tb_lattice coordinates.
        """
        for pos in self.__cluster_sites:
            yield pos

    def __repr__(self):
        def f(A):
            return list([ tuple(x) for x in A])
        return """SuperLattice class: \n
   Base TBLattice: %s
   SuperLattice Units: %s
   Remove internal Hoppings: %s
   Cluster site positions: %s"""%(self.__BaseLattice, f(self.__super_lattice_units), self.__cluster_sites, self.__remove_internal_hoppings)
