
import numpy as np

from pytriqs.archive import HDFArchive

from pytriqs.applications.susceptibility.Dummy import Dummy

d = Dummy()
d.test_data = np.random.random(100)

filename = 'Dummy.out.h5'

with HDFArchive(filename,'w') as ar:
    ar['d'] = d
    
with HDFArchive(filename,'r') as ar:
    d_ref = ar['d']

np.testing.assert_array_almost_equal(d.test_data, d_ref.test_data)
