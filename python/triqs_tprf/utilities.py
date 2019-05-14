# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 S. Käser
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

""" Collection of functions that can be useful in scripts """

# ----------------------------------------------------------------------

import subprocess
import os
import tarfile
from tempfile import NamedTemporaryFile

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

def show_version_info(info):
    """ Return a string that formats the version information

    Parameter:
    info : tuple of strings, coming from triqs_tprf.version.info
    """
    string = "TPRF version %s of hash %s and TRIQS version %s of hash %s"%info
    return string

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

def read_TarGZ_HDFArchive(filename):

    tar = tarfile.open(filename, "r:gz")
    f = tar.extractfile(tar.getmembers()[0])

    tmp = NamedTemporaryFile(delete=False)
    tmp.write(f.read())
    tmp.close()

    data = HDFArchive(tmp.name, 'r')

    os.remove(tmp.name)

    return data

def beta_to_temperature(beta):
    """Convert beta in 1/eV to Temperature in Kelvin
    """

    def eV_to_Kelvin(ev):
        return 11604.5250061657 * ev

    T = 1. / beta
    return eV_to_Kelvin(T)

def temperature_to_beta(T):
    """Convert Temperature in Kelvin to beta in 1/eV
    """

    def Kelvin_to_eV(K):
        return K / 11604.5250061657

    T = Kelvin_to_eV(T)
    beta = 1./ T
    return beta
