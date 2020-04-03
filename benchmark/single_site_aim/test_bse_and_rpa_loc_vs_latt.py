# ----------------------------------------------------------------------

""" Local and lattice RPA and BSE calculations comparison for
single site SIAM.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import os
import numpy as np

from pytriqs.gf import Gf, Idx
from pytriqs.operators import c, c_dag
from h5 import HDFArchive

from pytriqs.gf import MeshBrillouinZone, MeshProduct, MeshImFreq
from pytriqs.lattice.lattice_tools import BravaisLattice, BrillouinZone

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.freq_conv import block_iw_AB_to_matrix_valued
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH, chi_from_gg2_PH
from triqs_tprf.linalg import inverse_PH

# ----------------------------------------------------------------------
def trace_nn(G2_wnn):
    bmesh = G2_wnn.mesh.components[0]
    G2_w = Gf(mesh=bmesh, target_shape=G2_wnn.target_shape)
    G2_w.data[:] = np.sum(G2_wnn.data, axis=(1, 2)) / bmesh.beta**2
    return G2_w
    
# ----------------------------------------------------------------------
def read_TarGZ_HDFArchive(filename):

    import tarfile
    from tempfile import NamedTemporaryFile

    tar = tarfile.open(filename, "r:gz")
    f = tar.extractfile(tar.getmembers()[0])

    tmp = NamedTemporaryFile(delete=False)
    tmp.write(f.read())
    tmp.close()

    with HDFArchive(tmp.name, 'r') as res:
        p = res['p']

    os.remove(tmp.name)

    return p
    
# ----------------------------------------------------------------------
def make_calc():
            
    # ------------------------------------------------------------------
    # -- Read precomputed ED data

    filename = "data_pomerol.tar.gz"
    p = read_TarGZ_HDFArchive(filename)
    
    # ------------------------------------------------------------------
    # -- RPA tensor
    
    from triqs_tprf.rpa_tensor import get_rpa_tensor
    from triqs_tprf.rpa_tensor import fundamental_operators_from_gf_struct
    
    fundamental_operators = fundamental_operators_from_gf_struct(p.gf_struct)
    p.U_abcd = get_rpa_tensor(p.H_int, fundamental_operators)

    # ------------------------------------------------------------------
    # -- Generalized PH susceptibility
            
    loc_bse = ParameterCollection()
         
    loc_bse.chi_wnn = chi_from_gg2_PH(p.G_iw, p.G2_iw_ph)
    loc_bse.chi0_wnn = chi0_from_gg2_PH(p.G_iw, p.G2_iw_ph)
    
    loc_bse.gamma_wnn = inverse_PH(loc_bse.chi0_wnn) - inverse_PH(loc_bse.chi_wnn)
    loc_bse.chi_wnn_ref = inverse_PH( inverse_PH(loc_bse.chi0_wnn) - loc_bse.gamma_wnn )

    np.testing.assert_array_almost_equal(
        loc_bse.chi_wnn.data, loc_bse.chi_wnn_ref.data)
    
    loc_bse.chi0_w = trace_nn(loc_bse.chi0_wnn)
    loc_bse.chi_w = trace_nn(loc_bse.chi_wnn)

    # ------------------------------------------------------------------
    # -- RPA, using BSE inverses and constant Gamma

    loc_rpa = ParameterCollection()
    loc_rpa.U_abcd = p.U_abcd
    
    # -- Build constant gamma
    loc_rpa.gamma_wnn = loc_bse.gamma_wnn.copy()
    loc_rpa.gamma_wnn.data[:] = loc_rpa.U_abcd[None, None, None, ...]
    # Nb! In the three frequency form $\Gamma \propto U/\beta^2$
    loc_rpa.gamma_wnn.data[:] /= p.beta**2 

    loc_rpa.chi0_wnn = loc_bse.chi0_wnn
    loc_rpa.chi0_w = loc_bse.chi0_w
    
    # -- Solve RPA
    loc_rpa.chi_wnn = inverse_PH( inverse_PH(loc_rpa.chi0_wnn) - loc_rpa.gamma_wnn )
    loc_rpa.chi_w = trace_nn(loc_rpa.chi_wnn)
    
    # ------------------------------------------------------------------
    # -- Bubble RPA on lattice

    lat_rpa = ParameterCollection()
    
    # -- Setup dummy lattice Green's function equal to local Green's function
    
    bz = BrillouinZone(BravaisLattice(units=np.eye(3), orbital_positions=[(0,0,0)]))
    periodization_matrix = np.diag(np.array(list([1]*3), dtype=np.int32))
    kmesh = MeshBrillouinZone(bz, periodization_matrix)    
    wmesh = MeshImFreq(beta=p.beta, S='Fermion', n_max=p.nwf_gf)

    lat_rpa.g_wk = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=p.G_iw.target_shape)
    lat_rpa.g_wk[:, Idx(0, 0, 0)] = p.G_iw

    # -- chi0_wk bubble and chi_wk_rpa bubble RPA

    from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
    lat_rpa.chi0_wk = imtime_bubble_chi0_wk(lat_rpa.g_wk, nw=1)

    from triqs_tprf.lattice import solve_rpa_PH
    lat_rpa.chi_wk = solve_rpa_PH(lat_rpa.chi0_wk, p.U_abcd)

    lat_rpa.chi0_w = lat_rpa.chi0_wk[:, Idx(0,0,0)]
    lat_rpa.chi_w = lat_rpa.chi_wk[:, Idx(0,0,0)]

    print '--> cf Tr[chi0] and chi0_wk'
    print loc_rpa.chi0_w.data.reshape((4, 4)).real
    print lat_rpa.chi0_w.data.reshape((4, 4)).real

    np.testing.assert_array_almost_equal(
        loc_rpa.chi0_w.data, lat_rpa.chi0_w.data, decimal=2)

    print 'ok!'

    print '--> cf Tr[chi_rpa] and chi_wk_rpa'
    print loc_rpa.chi_w.data.reshape((4, 4)).real
    print lat_rpa.chi_w.data.reshape((4, 4)).real

    np.testing.assert_array_almost_equal(
        loc_rpa.chi_w.data, lat_rpa.chi_w.data, decimal=2)

    print 'ok!'
    
    # ------------------------------------------------------------------
    # -- Lattice BSE

    lat_bse = ParameterCollection()

    lat_bse.g_wk = lat_rpa.g_wk
    
    from triqs_tprf.lattice import fourier_wk_to_wr
    lat_bse.g_wr = fourier_wk_to_wr(lat_bse.g_wk)

    from triqs_tprf.lattice import chi0r_from_gr_PH
    lat_bse.chi0_wnr = chi0r_from_gr_PH(nw=1, nnu=p.nwf, gr=lat_bse.g_wr)

    from triqs_tprf.lattice import chi0q_from_chi0r
    lat_bse.chi0_wnk = chi0q_from_chi0r(lat_bse.chi0_wnr)

    # -- Lattice BSE calc
    from triqs_tprf.lattice import chiq_from_chi0q_and_gamma_PH
    lat_bse.chi_kwnn = chiq_from_chi0q_and_gamma_PH(lat_bse.chi0_wnk, loc_bse.gamma_wnn)
    
    # -- Trace results
    from triqs_tprf.lattice import chi0q_sum_nu_tail_corr_PH
    from triqs_tprf.lattice import chi0q_sum_nu
    lat_bse.chi0_wk_tail_corr = chi0q_sum_nu_tail_corr_PH(lat_bse.chi0_wnk)
    lat_bse.chi0_wk = chi0q_sum_nu(lat_bse.chi0_wnk)

    from triqs_tprf.lattice import chiq_sum_nu, chiq_sum_nu_q
    lat_bse.chi_kw = chiq_sum_nu(lat_bse.chi_kwnn)

    lat_bse.chi0_w_tail_corr = lat_bse.chi0_wk_tail_corr[:, Idx(0, 0, 0)]
    lat_bse.chi0_w = lat_bse.chi0_wk[:, Idx(0, 0, 0)]
    lat_bse.chi_w = lat_bse.chi_kw[Idx(0, 0, 0), :]
    
    print '--> cf Tr[chi0_wnk] and chi0_wk'
    print lat_bse.chi0_w_tail_corr.data.reshape((4, 4)).real
    print lat_bse.chi0_w.data.reshape((4, 4)).real
    print lat_rpa.chi0_w.data.reshape((4, 4)).real

    np.testing.assert_array_almost_equal(
        lat_bse.chi0_w_tail_corr.data, lat_rpa.chi0_w.data)

    np.testing.assert_array_almost_equal(
        lat_bse.chi0_w.data, lat_rpa.chi0_w.data, decimal=2)
    
    print 'ok!'
    
    print '--> cf Tr[chi_kwnn] and chi_wk'
    print lat_bse.chi_w.data.reshape((4, 4)).real
    print loc_bse.chi_w.data.reshape((4, 4)).real

    np.testing.assert_array_almost_equal(
        lat_bse.chi_w.data, loc_bse.chi_w.data)

    print 'ok!'
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_bse_rpa.h5'
    with HDFArchive(filename,'w') as res:
        res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc()
