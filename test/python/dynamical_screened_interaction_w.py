
# ----------------------------------------------------------------------

import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice

from triqs_tprf.lattice import lindhard_chi00


from triqs_tprf.gw import dynamical_screened_interaction_W

from triqs.gf import Gf, MeshImFreq, MeshReFreq, inverse, Idx
from triqs.gf.meshes import MeshDLRImFreq
from triqs.gf.mesh_product import MeshProduct

# ----------------------------------------------------------------------

def test_dynamical_screening_functions_single_orbital():
    print("Single orbital")
    nw = 100
    wmin = -5.0
    wmax = +5.0
    norb = 1
    beta = 1000.0
    mu = 2.0
    delta = 0.01
    
    t = -1.0 * np.eye(norb)
    
    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            (+1,) : t,
            (-1,) : t,
            },
        orbital_positions = [(0,0,0)]*norb,
        )
    
    kmesh = t_r.get_kmesh(n_k=(8, 1, 1))
    e_k = t_r.fourier(kmesh)
    
    
    kmesh = e_k.mesh
    wmesh = MeshImFreq(beta, 'Boson', nw)
    fmesh = MeshReFreq(wmin, wmax, nw)
    
    V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_k.data[:,0,0,0,0] = 1.0 / (np.linalg.norm(list(kmesh.values()), axis=1) + 0.01)
    
    V_wk = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=[norb]*4)
    for w in wmesh:
        wi = w.data_index
        V_wk.data[wi, :] = V_k.data[:]
    
    V_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=[norb]*4)
    for f in fmesh:
        fi = f.data_index
        V_fk.data[fi, :] = V_k.data[:]
    
    
    print('--> pi_bubble')
    PI_wk = lindhard_chi00(e_k=e_k, mesh=wmesh, mu=mu)
    PI_fk = lindhard_chi00(e_k=e_k, mesh=fmesh, beta=beta, mu=mu, delta=delta)

    print('--> reference screened_interaction_W')
    Wr_wk_ref = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=V_k.target_shape)
    for w in wmesh:
        wi = w.data_index
        Wr_wk_ref.data[wi,:] = V_k.data[:] / (1.0 - PI_wk.data[wi,:] * V_k.data[:])
    
    Wr_fk_ref = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=V_k.target_shape)
    for f in fmesh:
        fi = f.data_index
        Wr_fk_ref.data[fi,:] = V_k.data[:] / (1.0 - PI_fk.data[fi,:] * V_k.data[:])
    
    
    print('--> screened_interaction_W (Matsubara, static)')
    Wr_wk_1 = dynamical_screened_interaction_W(PI_wk, V_k)
    np.testing.assert_array_almost_equal(Wr_wk_1.data, Wr_wk_ref.data)
    
    print('--> screened_interaction_W (Matsubara, dynamic)')
    Wr_wk_2 = dynamical_screened_interaction_W(PI_wk, V_wk)
    np.testing.assert_array_almost_equal(Wr_wk_2.data, Wr_wk_ref.data)
    
    
    print('--> screened_interaction_W (real freq., static)')
    Wr_fk_1 = dynamical_screened_interaction_W(PI_fk, V_k)
    np.testing.assert_array_almost_equal(Wr_fk_1.data[:], Wr_fk_ref.data[:])
    
    print('--> screened_interaction_W (real freq., dynamic)')
    Wr_fk_2 = dynamical_screened_interaction_W(PI_fk, V_fk)
    np.testing.assert_array_almost_equal(Wr_fk_2.data[:], Wr_fk_ref.data[:])

def test_dynamical_screening_functions_multiple_orbitals():
    print("2 orbtials")
    nw = 100
    wmin = -5.0
    wmax = +5.0
    norb = 2
    beta = 1000.0
    mu = 2.0
    delta = 0.01
    
    lamb = 10.
    eps = 1e-8

    t = -1.0 * np.eye(norb)
    
    t_r = TBLattice(
        units = [(1, 0, 0)],
        hopping = {
            (+1,) : t,
            (-1,) : t,
            },
        orbital_positions = [(0,0,0)]*norb,
        )
    
    kmesh = t_r.get_kmesh(n_k=(8, 1, 1))
    e_k = t_r.fourier(kmesh)
    
    
    kmesh = e_k.mesh
    wmesh = MeshImFreq(beta, 'Boson', nw)
    DLRwmesh = MeshDLRImFreq(beta, 'Boson', lamb, eps)
    fmesh = MeshReFreq(wmin, wmax, nw)
    
    V_k = Gf(mesh=kmesh, target_shape=[norb]*4)
    V_k.data[:,0,0,0,0] = 5.0 / (np.linalg.norm(list(kmesh.values()), axis=1) + 0.01)
    V_k.data[:,0,1,0,1] = 2.0 / (np.linalg.norm(list(kmesh.values()), axis=1) + 0.05)
    V_k.data[:,1,0,1,0] = 1.0 / (np.linalg.norm(list(kmesh.values()), axis=1) + 0.1)
    V_k.data[:,1,1,1,1] = 7.0 / (np.linalg.norm(list(kmesh.values()), axis=1) + 0.06)
    
    V_wk = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=[norb]*4)
    for w in wmesh:
        wi = w.data_index
        V_wk.data[wi, :] = V_k.data[:]
    
    V_Dwk = Gf(mesh=MeshProduct(DLRwmesh, kmesh), target_shape=[norb]*4)
    for w in DLRwmesh:
        wi = w.data_index
        V_Dwk.data[wi, :] = V_k.data[:]

    V_fk = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=[norb]*4)
    for f in fmesh:
        fi = f.data_index
        V_fk.data[fi, :] = V_k.data[:]
    
    
    print('--> pi_bubble')
    PI_wk = lindhard_chi00(e_k=e_k, mesh=wmesh, mu=mu)
    PI_wk.data[:] = 0.0
    PI_wk.data[:,:,0,0,0,0] = -5.0
    PI_wk.data[:,:,1,0,1,0] = -2.0
    PI_wk.data[:,:,0,1,0,1] = -23.0
    PI_wk.data[:,:,1,1,1,1] = -12.0
    
    PI_Dwk = Gf(mesh=MeshProduct(DLRwmesh, kmesh), target_shape=PI_wk.target_shape)
    PI_Dwk.data[:] = 0.0
    PI_Dwk.data[:,:,0,0,0,0] = -5.0
    PI_Dwk.data[:,:,1,0,1,0] = -2.0
    PI_Dwk.data[:,:,0,1,0,1] = -23.0
    PI_Dwk.data[:,:,1,1,1,1] = -12.0

    PI_fk = lindhard_chi00(e_k=e_k, mesh=fmesh, beta=beta, mu=mu, delta=delta)
    PI_fk.data[:] = 0.0
    PI_fk.data[:,:,0,0,0,0] = -5.0
    PI_fk.data[:,:,1,0,1,0] = -2.0
    PI_fk.data[:,:,0,1,0,1] = -23.0
    PI_fk.data[:,:,1,1,1,1] = -12.0
    
    print('--> reference screened_interaction_W')
    Wr_wk_ref = Gf(mesh=MeshProduct(wmesh, kmesh), target_shape=V_k.target_shape)
    for w in wmesh:
        wi = w.data_index
        Wr_wk_ref.data[wi,:] = V_k.data[:] / (1.0 - PI_wk.data[wi,:] * V_k.data[:])
    
    Wr_Dwk_ref = Gf(mesh=MeshProduct(DLRwmesh, kmesh), target_shape=V_k.target_shape)
    for w in DLRwmesh:
        wi = w.data_index
        Wr_Dwk_ref.data[wi,:] = V_k.data[:] / (1.0 - PI_Dwk.data[wi,:] * V_k.data[:])

    Wr_fk_ref = Gf(mesh=MeshProduct(fmesh, kmesh), target_shape=V_k.target_shape)
    for f in fmesh:
        fi = f.data_index
        Wr_fk_ref.data[fi,:] = V_k.data[:] / (1.0 - PI_fk.data[fi,:] * V_k.data[:])
    
    
    print('--> screened_interaction_w (matsubara, static)')
    Wr_wk_1 = dynamical_screened_interaction_W(PI_wk, V_k)
    np.testing.assert_array_almost_equal(Wr_wk_1.data, Wr_wk_ref.data)
    
    print('--> screened_interaction_W (Matsubara, dynamic)')
    Wr_wk_2 = dynamical_screened_interaction_W(PI_wk, V_wk)
    np.testing.assert_array_almost_equal(Wr_wk_2.data, Wr_wk_ref.data)
    
    print('--> screened_interaction_w (DLR matsubara, static)')
    Wr_Dwk_1 = dynamical_screened_interaction_W(PI_Dwk, V_k)
    np.testing.assert_array_almost_equal(Wr_Dwk_1.data, Wr_Dwk_ref.data)

    print('--> screened_interaction_W (DLR Matsubara, dynamic)')
    Wr_Dwk_2 = dynamical_screened_interaction_W(PI_Dwk, V_Dwk)
    np.testing.assert_array_almost_equal(Wr_Dwk_2.data, Wr_Dwk_ref.data)
    
    print('--> screened_interaction_W (real freq., static)')
    Wr_fk_1 = dynamical_screened_interaction_W(PI_fk, V_k)
    np.testing.assert_array_almost_equal(Wr_fk_1.data[:], Wr_fk_ref.data[:])
    
    print('--> screened_interaction_W (real freq., dynamic)')
    Wr_fk_2 = dynamical_screened_interaction_W(PI_fk, V_fk)
    np.testing.assert_array_almost_equal(Wr_fk_2.data[:], Wr_fk_ref.data[:])

 
    
    
if __name__ == "__main__":
    test_dynamical_screening_functions_single_orbital()
    test_dynamical_screening_functions_multiple_orbitals() 
    
