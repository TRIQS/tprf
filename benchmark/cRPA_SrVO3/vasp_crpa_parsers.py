"""

Author: H. U.R. Strand
"""

import glob
import warnings
import numpy as np

def read_vasp_crpa_momentum_space_interaction(path, prefix, verbose=False):

    filenames = glob.glob(path + '/' + prefix + '*')
    filenames = np.sort(filenames)

    data = np.zeros((0, 9))
    for filename in filenames:
        if verbose: print '--> Loading:', filename
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

    if verbose: print 'norb, nq =', norb, nq

    Q = np.zeros((nq, 3), dtype=np.float)
    U_Q = np.zeros((nq, norb, norb, norb, norb), dtype=np.complex)
    
    for qidx in xrange(nq):
        u = np.zeros([norb]*4, dtype=np.complex)
        s, e = qidx * norb**4, (qidx + 1) * norb**4
        np.testing.assert_array_almost_equal(q[s:e] - q[s], np.zeros_like(q[s:e]))
        ijkl, v = idxs[s:e], vals[s:e]
        for n in xrange(norb**4):
            u[tuple(ijkl[n] - 1)] = v[n]

        Q[qidx] = q[s]
        U_Q[qidx] = u

    Q, uidx = np.unique(Q, return_index=True, axis=0)
    U_Q = U_Q[uidx]
        
    return U_Q, Q
