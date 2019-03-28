# ----------------------------------------------------------------------

""" Collection of functions that can be useful for testing scripts """

# ----------------------------------------------------------------------

import subprocess
import os
import tarfile
from tempfile import NamedTemporaryFile

# ----------------------------------------------------------------------

from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

def get_git_revision_short_hash():
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip()

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
