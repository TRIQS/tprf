# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 S. Käser
# Copyright (C) 2019 by The Simons Foundation
# Author: H. U.R. Strand
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

# ----------------------------------------------------------------------
def show_version_info(info):
    """ Return a string that formats the version information

    Parameter:
    info : tuple of strings, coming from triqs_tprf.version.info
    """
    string = "TPRF version %s of hash %s and TRIQS version %s of hash %s"%info
    return string

# ----------------------------------------------------------------------
def write_TarGZ_HDFArchive(filename, **kwargs):

    import os
    import tarfile
    from h5 import HDFArchive
    
    filename = filename.split('.')[0]
    filename_h5 = filename + '.h5'
    filename_tar = filename + '.tar.gz'

    with HDFArchive(filename_h5, 'w') as res:
        for key, value in kwargs.items():
            res[key] = value

    with tarfile.open(filename_tar, 'w:gz') as tar:
        tar.add(filename_h5)

    os.remove(filename_h5)

# ----------------------------------------------------------------------
def read_TarGZ_HDFArchive(filename):

    import os
    import tarfile
    from tempfile import NamedTemporaryFile
    from h5 import HDFArchive

    tar = tarfile.open(filename, "r:gz")
    f = tar.extractfile(tar.getmembers()[0])

    tmp = NamedTemporaryFile(delete=False)
    tmp.write(f.read())
    tmp.close()

    data = HDFArchive(tmp.name, 'r')

    os.remove(tmp.name)

    return data

# ----------------------------------------------------------------------
def BlockGf_data(G):
    """ Returns a ndarray copy of all data in a BlockGf """

    import numpy as np
    
    shape = [G.n_blocks] + list(G[G.indices.next()].data.shape)
    data = np.zeros(shape, dtype=np.complex)
    for bidx, (b, g) in enumerate(G):
        data[bidx] = g.data.copy()

    return data

# ----------------------------------------------------------------------
def legendre_filter(G_tau, order=100, G_l_cut=1e-19):
    """ Filter binned imaginary time Green's function
    using a Legendre filter of given order and coefficient threshold. 
    
    Parameters
    ----------

    G_tau : TRIQS imaginary time Block Green's function

    order : int
        Legendre expansion order in the filter

    G_l_cut : float
        Legendre coefficient cut-off 

    Returns
    -------

    G_l : TRIQS Legendre Block Green's function
        Fitted Green's function on a Legendre mesh

    """
    
    import numpy as np
    from pytriqs.gf import BlockGf
    from pytriqs.gf.tools import fit_legendre
    from pytriqs.gf.gf_fnt import enforce_discontinuity
    
    l_g_l = []

    for b, g in G_tau:

        g_l = fit_legendre(g, order=order)
        g_l.data[:] *= (np.abs(g_l.data) > G_l_cut)
        enforce_discontinuity(g_l, np.array([[1.]]))
        l_g_l.append(g_l)

    G_l = BlockGf(name_list=list(G_tau.indices), block_list=l_g_l)
    return G_l

# ----------------------------------------------------------------------
def G2_loc_fixed_fermionic_window_python(g2, nwf):

    """ Limit the last two fermionic freqiencies of a three
    frequency Green's function object :math:`G(\omega, \nu, \nu')`
    to ``nwf``. """

    nw = (g2.data.shape[0] + 1) / 2
    n = g2.data.shape[1]
    beta = g2.mesh.components[0].beta
    
    assert(n/2 >= nwf)

    from pytriqs.gf import Gf, MeshImFreq, MeshProduct
    
    mesh_iw = MeshImFreq(beta=beta, S='Boson', n_max=nw)
    mesh_inu = MeshImFreq(beta=beta, S='Fermion', n_max=nwf)
    mesh_prod = MeshProduct(mesh_iw, mesh_inu, mesh_inu)

    g2_out = Gf(mesh=mesh_prod, target_shape=g2.target_shape)

    s = n/2 - nwf
    e = n/2 + nwf
    
    g2_out.data[:] = g2.data[:, s:e, s:e]
    return g2_out
