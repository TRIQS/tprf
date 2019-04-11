#!/usr/bin/python

import numpy as np
import collections
import itertools

'''
python functions for reading Uijkl from a VASP cRPA run and the evaluating the matrix
elements for different basis sets.

Copyright (C) 2019, A. Hampel and C. Ederer from Materials Theory Group
at ETH Zurich

Authors: A. Hampel, C. Ederer, H. U.R. Strand
'''

def read_uijkl(path_to_uijkl):
    '''
    reads the VASP UIJKL files or the vijkl file if wanted

    Parameters
    ----------
    path_to_uijkl : string
        path to Uijkl like file
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom

    __Returns:__
    uijkl : numpy array
        uijkl Coulomb tensor

    '''
    data = np.loadtxt(path_to_uijkl)
    dim = int(data[:, 0].max())
    uijkl = np.zeros((dim,dim,dim,dim), dtype=np.complex)
 
    for line in range(data.shape[0]):
        idxs = tuple(np.array(data[line, :4], dtype=np.int) - 1)
        uijkl[idxs] = data[line,4] + 1.j * data[line, 5]

    return uijkl

def red_to_2ind(uijkl, verbose=False):
    '''
    reduces the 4index coulomb matrix to a 2index matrix and
    follows the procedure given in PRB96 seth,peil,georges:
    m = ii , m'=jj
    U_antipar = U_mm'^oo' = U_mm'mm' (Coulomb Int)
    U_par = U_mm'^oo = U_mm'mm' - U_mm'm'm (for intersite interaction)
    U_ijij (Hunds coupling)

    the indices in VASP are switched: U_ijkl ---VASP--> U_ikjl 

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    __Returns:__
    Uij_anti : numpy array
        red 2 index matrix U_mm'mm'
    Uijij : numpy array
        red 2 index matrix U_ijij (Hunds coupling)
    Uijji : numpy array
        red 2 index matrix Uijji
    Uij_par : numpy array
        red 2 index matrix U_mm\'mm\' - U_mm\'m\'m

    '''
    assert( (np.array(uijkl.shape) - uijkl.shape[0] == np.zeros(4)).all() )
    dim = uijkl.shape[0]

    # create 2 index matrix
    Uij_anti = np.zeros((dim, dim), dtype=uijkl.dtype)
    Uij_par = np.zeros_like(Uij_anti)
    Uijij = np.zeros_like(Uij_anti)
    Uijji = np.zeros_like(Uij_anti)

    for i, j in itertools.product(range(dim), repeat=2):
        # the indices in VASP are switched: U_ijkl ---VASP--> U_ikjl
        Uij_anti[i, j] = uijkl[i, i, j, j]
        Uijij[i, j] = uijkl[i, j, i, j]
        Uijji[i, j] = uijkl[i, j, j, i]
        Uij_par[i, j] = uijkl[i, i, j, j] - uijkl[i, j, j, i]

    if verbose:
        np.set_printoptions(precision=3,suppress=True)
        print 'reduced U anti-parallel = U_mm\'\^oo\' = U_mm\'mm\' matrix : \n', Uij_anti
        print 'reduced Uijij : \n', Uijij
        print 'reduced Uijji : \n', Uijji
        print 'reduced U parallel = U_mm\'\^oo = U_mm\'mm\' - U_mm\'m\'m matrix : \n', Uij_par

    return Uij_anti, Uijij, Uijji, Uij_par


def calc_kan_params_egeg(uijkl,n_sites,n_orb,out=False):
    '''
    calculates the kanamori interaction parameters from a
    given Uijkl matrix. Follows the procedure given in
    PHYSICAL REVIEW B 86, 165105 (2012) Vaugier,Biermann
    formula 30,31,32

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    __Returns:__
    int_params : direct
        kanamori parameters
    '''

    int_params = collections.OrderedDict()
    dim = n_sites*n_orb
    Uij_anti,Uijij,Uij_par = red_to_2ind(uijkl,n_sites,n_orb,out)

    # calculate intra-orbital U
    U = 0.0
    for i in range(0,n_orb):
        U += uijkl[i,i,i,i]
    U = U/(n_orb)
    int_params['U'] = U

    # calculate the U'
    Uprime = 0.0
    for i in range(0,n_orb):
        for j in range(0,n_orb):
            if i != j:
                Uprime +=  uijkl[i,i,j,j]
    Uprime = Uprime/ (n_orb*(n_orb-1))
    int_params['Uprime'] = Uprime

    # calculate J
    J = 0.0
    for i in range(0,n_orb):
        for j in range(0,n_orb):
            if i != j:
                J +=  uijkl[i,j,i,j]
    J = J/ (n_orb*(n_orb-1))
    int_params['J'] = J

    if out:
        print 'U= ', "{:.4f}".format(U)
        print 'U\'= ', "{:.4f}".format(Uprime)
        print 'J= ', "{:.4f}".format(J)

    return int_params

def calc_u_avg_fulld(uijkl,n_sites,n_orb,out=False):
    '''
    calculates the coulomb integrals from a
    given Uijkl matrix for full d shells. Follows the procedure given
    in Pavarini  - 2014 - arXiv - 1411 6906 - julich school U matrix
    page 8 or as done in
    PHYSICAL REVIEW B 86, 165105 (2012) Vaugier,Biermann
    formula 23, 25

    works atm only for full d shell (l=2)

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    __Returns:__
    int_params : direct
        Slater parameters
    '''

    int_params = collections.OrderedDict()
    dim = n_sites*n_orb
    Uij_anti,Uijij,Uijji,Uij_par = red_to_2ind(uijkl,n_sites,n_orb,out=False)
    # U_antipar = U_mm'^oo' = U_mm'mm' (Coulomb Int)
    # U_par = U_mm'^oo = U_mm'mm' - U_mm'm'm (for intersite interaction)
    # U_ijij (Hunds coupling)


    # calculate intra-orbital U
    U = 0.0
    for i in range(0,n_orb):
        for j in range(0,n_orb):
            U += Uij_anti[i,j]
    # 1/(2l+1)^2
    U = U/(5*5)
    int_params['U'] = U

    # calculate J
    J = 0.0
    for i in range(0,n_orb):
        for j in range(0,n_orb):
            if i != j:
                J +=  Uijji[i,j]
    J = J/ (20)
    # 20 for 2l(2l+1)
    int_params['J'] = J

    if out:
        print 'U= ', "{:.4f}".format(U)
        print 'J= ', "{:.4f}".format(J)

    return int_params
