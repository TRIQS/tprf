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

import os
import tarfile
from tempfile import NamedTemporaryFile

import numpy as np

from pytriqs.archive import HDFArchive

from pytriqs.gf import Gf, MeshImFreq, MeshProduct, BlockGf
from pytriqs.gf.tools import fit_legendre
from pytriqs.gf.gf_fnt import enforce_discontinuity

from triqs_tprf.lattice import lattice_dyson_g0_wk, solve_rpa_PH, gamma_PP_singlet
from triqs_tprf.tight_binding import create_model_for_tests
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
from triqs_tprf.rpa_tensor import kanamori_charge_and_spin_quartic_interaction_tensors

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

    
    mesh_iw = MeshImFreq(beta=beta, S='Boson', n_max=nw)
    mesh_inu = MeshImFreq(beta=beta, S='Fermion', n_max=nwf)
    mesh_prod = MeshProduct(mesh_iw, mesh_inu, mesh_inu)

    g2_out = Gf(mesh=mesh_prod, target_shape=g2.target_shape)

    s = n/2 - nwf
    e = n/2 + nwf
    
    g2_out.data[:] = g2.data[:, s:e, s:e]
    return g2_out

# ----------------------------------------------------------------------
def beta_to_temperature(beta):
    """Convert beta in 1/eV to Temperature in Kelvin
    """
    def eV_to_Kelvin(ev):
        return 11604.5250061657 * ev

    T = 1. / beta
    return eV_to_Kelvin(T)

# ----------------------------------------------------------------------
def temperature_to_beta(T):
    """Convert Temperature in Kelvin to beta in 1/eV
    """
    def Kelvin_to_eV(K):
        return K / 11604.5250061657

    T = Kelvin_to_eV(T)
    beta = 1./ T
    return beta

# ----------------------------------------------------------------------
def create_eliashberg_ingredients(p):
    H = create_model_for_tests(**p)
    e_k = H.on_mesh_brillouin_zone(n_k=[p.nk] * p.dim + [1] * (3 - p.dim))

    wmesh = MeshImFreq(beta=p.beta, S="Fermion", n_max=p.nw)
    g0_wk = lattice_dyson_g0_wk(mu=p.mu, e_k=e_k, mesh=wmesh)

    chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw)

    U_c, U_s = kanamori_charge_and_spin_quartic_interaction_tensors(
        p.norb, p.U, p.Up, p.J, p.Jp
    )

    chi_s = solve_rpa_PH(chi0_wk, U_s)
    chi_c = solve_rpa_PH(chi0_wk, -U_c)  # Minus for correct charge rpa equation

    gamma = gamma_PP_singlet(chi_c, chi_s, U_c, U_s)

    eliashberg_ingredients = ParameterCollection(
                                g0_wk = g0_wk,
                                gamma = gamma,
                                U_s = U_s,
                                U_c = U_c,
                                chi_s = chi_s,
                                chi_c = chi_c,
                                )
    return eliashberg_ingredients

# ----------------------------------------------------------------------
def assert_parameter_collection_not_equal_model_parameters(p1, p2, model_parameters):
    for model_parameter in model_parameters:
        value1, value2 = p1[model_parameter], p2[model_parameter]
        if value1 != value2:
            error = 'The model of the benchmark and the one used now are not the same.\n' 
            error += '\t\tNow: {0} = {1}, benchmark: {0} = {2}.'.format(model_parameter, value1, value2)
            raise AssertionError, error
