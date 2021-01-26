# ----------------------------------------------------------------------

""" Test if calculating chi0_wnk only a a specififc w point gives the same
result as the function used for a whole bosonic Matsubara mesh. """


# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, Idx, MeshProduct
from triqs_tprf.utilities import create_g0_wk_for_test_model
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.bse import get_chi0_wnk, get_chi0_nk_at_specific_w

# ----------------------------------------------------------------------


def test_chi0_nk_at_specific_w_against_full(g0_wk, p):
    chi0_wnk = get_chi0_wnk(g0_wk, nw=p.nw_chi, nwf=p.nwf)
    chi0_nk_at_specific_w = get_chi0_nk_at_specific_w(
        g0_wk, nw_index=p.nw_index, nwf=p.nwf
    )

    np.testing.assert_allclose(
        chi0_wnk[Idx(p.nw_index), :, :].data, chi0_nk_at_specific_w.data
    )


def test_chi0_nk_at_specific_w_mesh_types(g0_wk, p):
    chi0_nk_at_specific_w = get_chi0_nk_at_specific_w(
        g0_wk, nw_index=p.nw_index, nwf=p.nwf
    )

    assert isinstance(chi0_nk_at_specific_w.mesh, MeshProduct)
    assert isinstance(chi0_nk_at_specific_w.mesh[0], MeshImFreq)


if __name__ == "__main__":
    p = ParameterCollection(
        dim=2,
        norb=2,
        t1=1.0,
        t2=0.5,
        t12=0.1,
        t21=0.1,
        mu=0.0,
        beta=1,
        nk=2,
        nw=20,
        nw_chi=10,
        nwf=10,
        nw_index=3,
    )

    g0_wk = create_g0_wk_for_test_model(p)

    test_chi0_nk_at_specific_w_against_full(g0_wk, p)
    test_chi0_nk_at_specific_w_mesh_types(g0_wk, p)
