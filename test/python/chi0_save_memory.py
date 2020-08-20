# ----------------------------------------------------------------------

""" Test if memory saving version of the bare bubble chi0 is implemented
correctly. """

# ----------------------------------------------------------------------

import itertools
import time

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.utilities import create_g0_wk_for_test_model
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------

def test_chi0_wk_save_memory(g0_wk, p):
    chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=p.nw_chi, save_memory=False)
    chi0_wk_save_memory = imtime_bubble_chi0_wk(g0_wk, nw=p.nw_chi, save_memory=True)

    assert np.allclose(chi0_wk.data, chi0_wk_save_memory.data)


if __name__ == "__main__":
    p = ParameterCollection(
        dim=2,
        norb=2,
        t1 = 1.0,
        t2 = 0.5,
        t12 = 0.1,
        t21 = 0.1,
        mu=0.0,
        beta=1,
        nk=16,
        nw=100,
        nw_chi=10,
    )

    g0_wk = create_g0_wk_for_test_model(p)
    
    test_chi0_wk_save_memory(g0_wk, p)
