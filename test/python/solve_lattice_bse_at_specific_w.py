# ----------------------------------------------------------------------

""" Test if calculating the lattice susceptibility via the Bethe-
    Salpeter equation for a specific \omega gives the same result as the
    function used for a whole bosonic Matsubara mesh. 

Author: Stefan Käser (2020) stefan.kaeser7@gmail.com """


# ----------------------------------------------------------------------

import itertools

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

from triqs.gf import Gf, MeshImFreq, Idx, MeshProduct, MeshBrZone
from triqs_tprf.utilities import create_g0_wk_for_test_model
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.bse import solve_lattice_bse, solve_lattice_bse_at_specific_w

# ----------------------------------------------------------------------


def create_random_gamma_wnn(p):
    wmesh_gamma = MeshImFreq(beta=p.beta, S="Boson", n_max=p.nw_gamma)
    nmesh_gamma = MeshImFreq(beta=p.beta, S="Fermion", n_max=p.nwf)

    gamma_wnn = Gf(
        mesh=MeshProduct(wmesh_gamma, nmesh_gamma, nmesh_gamma),
        target_shape=2 * g0_wk.target_shape,
    )

    np.random.seed(p.seed)
    gamma_wnn.data[:] = np.random.rand(*gamma_wnn.data.shape)

    return gamma_wnn


def test_solve_lattice_bse_at_specific_w_against_full(g0_wk, gamma_wnn, nw_index):
    chi_kw, chi0_kw = solve_lattice_bse(g0_wk, gamma_wnn)
    chi_k_at_specific_w, chi0_k_at_specific_w = solve_lattice_bse_at_specific_w(
        g0_wk, gamma_wnn, nw_index=nw_index
    )

    np.testing.assert_allclose(
        chi0_kw[:, Idx(nw_index)].data, chi0_k_at_specific_w.data, atol=10e-16
    )
    np.testing.assert_allclose(
        chi_kw[:, Idx(nw_index)].data, chi_k_at_specific_w.data, atol=10e-16
    )


def test_lattice_bse_at_specific_w_mesh_types(g0_wk, gamma_wnn, nw_index):
    chi_k_at_specific_w, chi0_k_at_specific_w = solve_lattice_bse_at_specific_w(
        g0_wk, gamma_wnn, nw_index=nw_index
    )

    assert isinstance(chi_k_at_specific_w.mesh, MeshBrZone)
    assert isinstance(chi0_k_at_specific_w.mesh, MeshBrZone)


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
        nw=100,
        nw_gamma=10,
        nwf=10,
        nw_index=3,
        seed=101,
    )

    g0_wk = create_g0_wk_for_test_model(p)
    gamma_wnn = create_random_gamma_wnn(p)

    test_solve_lattice_bse_at_specific_w_against_full(g0_wk, gamma_wnn, p.nw_index)
    test_lattice_bse_at_specific_w_mesh_types(g0_wk, gamma_wnn, p.nw_index)
    
