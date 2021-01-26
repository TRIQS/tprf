# ----------------------------------------------------------------------

""" Test if accessing chi0_wnr only a a specififc w point gives the same
result as the function used for a whole bosonic Matsubara mesh. 

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """

# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import MeshImFreq, Idx, MeshProduct
from triqs_tprf.utilities import create_g0_wk_for_test_model
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi0r_from_gr_PH, chi0_nr_from_gr_PH_at_specific_w

# ----------------------------------------------------------------------


def test_chi0_nr_at_specific_w_against_full(g0_wr, p):
    chi0_wnr = chi0r_from_gr_PH(nw=p.nw_chi, nn=p.nwf, g_nr=g0_wr)
    chi0_nr_at_specific_w = chi0_nr_from_gr_PH_at_specific_w(
        nw_index=p.nw_index, nn=p.nwf, g_nr=g0_wr
    )

    np.testing.assert_allclose(
        chi0_wnr[Idx(p.nw_index), :, :].data, chi0_nr_at_specific_w.data
    )


def test_chi0_nr_at_specific_w_mesh_types(g0_wr, p):
    chi0_nr_at_specific_w = chi0_nr_from_gr_PH_at_specific_w(
        nw_index=p.nw_index, nn=p.nwf, g_nr=g0_wr
    )

    assert isinstance(chi0_nr_at_specific_w.mesh, MeshProduct)
    assert isinstance(chi0_nr_at_specific_w.mesh[0], MeshImFreq)


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
    g0_wr = fourier_wk_to_wr(g0_wk)

    test_chi0_nr_at_specific_w_against_full(g0_wr, p)
    test_chi0_nr_at_specific_w_mesh_types(g0_wr, p)
