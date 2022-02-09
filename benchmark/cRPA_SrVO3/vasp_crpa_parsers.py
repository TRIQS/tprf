"""

Author: H. U.R. Strand
"""

import glob
import warnings
import numpy as np

def read_vasp_crpa_momentum_space_interaction_to_ndarray(path, prefix, verbose=False):

    filenames = glob.glob(path + '/' + prefix + '*')
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
    nq = q.shape[0] / norb**4

    if verbose: print('norb, nq =', norb, nq)

    Q = np.zeros((nq, 3), dtype=np.float)
    U_Q = np.zeros((nq, norb, norb, norb, norb), dtype=complex)
    
    for qidx in range(nq):
        u = np.zeros([norb]*4, dtype=complex)
        s, e = qidx * norb**4, (qidx + 1) * norb**4
        np.testing.assert_array_almost_equal(q[s:e] - q[s], np.zeros_like(q[s:e]))
        ijkl, v = idxs[s:e], vals[s:e]
        for n in range(norb**4):
            u[tuple(ijkl[n] - 1)] = v[n]

        Q[qidx] = q[s]
        U_Q[qidx] = u

    Q, uidx = np.unique(Q, return_index=True, axis=0)
    U_Q = U_Q[uidx]
        
    return U_Q, Q

def convert_from_ndarray_to_triqs(U_Q, Q, cell, kpts):

    from triqs.gf import Gf, MeshBrZone
    from triqs.lattice.lattice_tools import BrillouinZone
    from triqs.lattice.lattice_tools import BravaisLattice

    bl = BravaisLattice(cell, [(0,0,0)])
    bz = BrillouinZone(bl)
    bzmesh = MeshBrZone(bz, np.diag(np.array(kpts, dtype=np.int32)))

    u_q = Gf(mesh=bzmesh, target_shape=U_Q.shape[1:])

    tmp = np.array(Q * kpts[None, :], dtype=np.int)
    I = [tuple(tmp[i]) for i in range(Q.shape[0])]

    for qidx, i in enumerate(I):
        for k in bzmesh:

            # -- Generalize this transform absolute to relative k-points
            a = cell[0, 0]
            q = k * 0.5 * a / np.pi

            j = tuple(np.array(kpts * q, dtype=np.int))
            if i == j: u_q[k].data[:] = U_Q[qidx]

    return u_q

def read_vasp_crpa_momentum_space_interaction_to_triqs(path, cell, kpts, verbose=False):

    UR_Q, Q = read_vasp_crpa_momentum_space_interaction_to_ndarray(path, 'UIJKL_Q_full.q*', verbose=verbose)
    VR_Q, Q = read_vasp_crpa_momentum_space_interaction_to_ndarray(path, 'VIJKL_Q_full.q*', verbose=verbose)
    VRR_Q, Q = read_vasp_crpa_momentum_space_interaction_to_ndarray(path, 'VIJKL_Q_redu.q*', verbose=verbose)
    U_Q = UR_Q + ( VR_Q - VRR_Q )
    
    u_q = convert_from_ndarray_to_triqs(U_Q, Q, cell, kpts)

    return u_q

