"""
Authors: H. U.R. Strand
         M. Roesner
"""

import glob
import warnings
import itertools
import numpy as np

def read_vasp_crpa_vq_to_ndarray(path, prefix, verbose=False, orbsub=[]):

    filenames = glob.glob(path + '/' + prefix)
    filenames = np.sort(filenames)

    data = np.zeros((0, 9))
    for filename in filenames:
        if verbose: print('--> Loading:', filename)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            d = np.genfromtxt(filename, skip_header=35)
        if len(d.shape) == 2:
            data = np.vstack((data, d))

    q = data[:, :3]
    idxs = np.array(data[:, 3:7], dtype=np.int)
    vals = data[:, 7] + 1.j * data[:, 8]

    norb = idxs.max()
    norbsub = len(orbsub)
    nq = int(q.shape[0] / norb**4)

    if verbose:
        print('  norb, nq =', norb, nq)
        if(norbsub > 0):
            print('orb subset =', orbsub)

    Q = np.zeros((nq, 3), dtype=np.float)

    if(norbsub == 0):
        norbsub = norb
        orbsub = range(norb)

    U_Q = np.zeros((nq, norbsub, norbsub, norbsub, norbsub), dtype=np.complex)
    
    for qidx in range(nq):
        u = np.zeros([norb]*4, dtype=np.complex)
        s, e = qidx * norb**4, (qidx + 1) * norb**4
        np.testing.assert_array_almost_equal(q[s:e] - q[s], np.zeros_like(q[s:e]))
        ijkl, v = idxs[s:e], vals[s:e]
        for n in range(norb**4):
            u[tuple(ijkl[n] - 1)] = v[n]

        for a, b, c, d in itertools.product(range(norbsub), repeat=4):
            # NOTE: See the different ordering of a,b,c,d on the left and
            #       right sides! This is due to the different odering
            #       between VASP and TPRF / TRIQS !!!!
            #
            #        U[a,b,c,d] --VASP--> U[b,a,d,c]
            #
            U_Q[qidx,a,b,c,d] = u[orbsub[b], orbsub[a], orbsub[d], orbsub[c]]

        Q[qidx] = q[s]

        if verbose: print('--> Reshaping q=[%+f, %+f, %+f] = %f %f'%(Q[qidx,0], Q[qidx,1], Q[qidx,2], np.real(u[0,0,0,0]), np.imag(u[0,0,0,0])))

    Q, uidx = np.unique(Q, return_index=True, axis=0)
    U_Q = U_Q[uidx]
        
    return U_Q, Q

def convert_from_ndarray_to_triqs(U_Q, Q, units, orbital_positions, kpts):

    from triqs.gf import Gf, MeshBrillouinZone
    from triqs.lattice.lattice_tools import BrillouinZone
    from triqs.lattice.lattice_tools import BravaisLattice
    from triqs_tprf.tight_binding import TBLattice
    from triqs_tprf.lattice_utils import get_relative_k_from_absolute

    bl = BravaisLattice(units, orbital_positions)
    bz = BrillouinZone(bl)
    bzmesh = MeshBrillouinZone(bz, np.diag(np.array(kpts, dtype=np.int32)))

    u_q = Gf(mesh=bzmesh, target_shape=U_Q.shape[1:])

    H  = TBLattice(units=units, hopping={}, orbital_positions=orbital_positions)
    kunits = H.bz.units()
   
    # note: np.rint() absolutely imporant to have correct
    #       rounding behaviour
    tmp = np.array(np.rint(Q * kpts[None, :]), dtype=np.int)
    I = [tuple(tmp[i]) for i in range(Q.shape[0])]

    for qidx, i in enumerate(I):

        i = np.array(i)
        
        if(i[0] < 0):
            i[0] += kpts[0]
        if(i[1] < 0):
            i[1] += kpts[1]
        if(i[2] < 0):
            i[2] += kpts[2]
            
        i = tuple(i)

        for k in bzmesh:

            q = get_relative_k_from_absolute(k.value, kunits)
            # note: np.rint() absolutely imporant to have correct
            #       rounding behaviour
            j = tuple(np.array(np.rint(kpts * q), dtype=np.int))
       
            if i == j: 
                u_q.data[k.linear_index,:] = U_Q[qidx]

    return u_q

def read_vasp_crpa_vq_to_triqs(path, units, orbital_positions, kpts, verbose=False, orbsub=[]):

    UR_Q, Q = read_vasp_crpa_vq_to_ndarray(path, 'UIJKL_Q_full.q*', verbose=verbose, orbsub=orbsub)
    VR_Q, Q = read_vasp_crpa_vq_to_ndarray(path, 'VIJKL_Q_full.q*', verbose=verbose, orbsub=orbsub)
    VRR_Q, Q = read_vasp_crpa_vq_to_ndarray(path, 'VIJKL_Q_redu.q*', verbose=verbose, orbsub=orbsub)
    U_Q = UR_Q + ( VR_Q - VRR_Q )
    
    v_q = convert_from_ndarray_to_triqs(VR_Q, Q, units, orbital_positions, kpts)
    u_q = convert_from_ndarray_to_triqs(U_Q, Q, units, orbital_positions, kpts)

    return v_q, u_q

