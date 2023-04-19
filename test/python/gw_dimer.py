
import itertools
import numpy as np
from scipy.linalg import block_diag

from triqs.lattice.tight_binding import TBLattice
from triqs.gf import Gf, MeshImFreq, Idx, MeshImTime
from triqs.operators import n, c, c_dag, Operator, dagger

#from triqs_tprf.rpa_tensor import get_rpa_tensor
from triqs_tprf.OperatorUtils import get_operator_index_map
from triqs_tprf.OperatorUtils import quadratic_matrix_from_operator
from triqs_tprf.OperatorUtils import operator_single_particle_transform

from triqs_tprf.hf_solver import HartreeSolver
from triqs_tprf.hf_solver import HartreeFockSolver

from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk
from triqs_tprf.lattice import split_into_dynamic_wk_and_constant_k

from triqs_tprf.gw import bubble_PI_wk
from triqs_tprf.gw import dynamical_screened_interaction_W
from triqs_tprf.lattice import gw_sigma, hartree_sigma, fock_sigma
from triqs_tprf.gw import g0w_sigma

from triqs.gf import Gf, MeshImFreq, Idx, MeshImTime
from triqs.gf.mesh_product import MeshProduct


def print_tensor(U, tol=1e-9):
    assert( len(U.shape) == 4)
    n = U.shape[0]
    
    for i,j,k,l in itertools.product(range(n), repeat=4):
        value = U[i, j, k, l]
        if np.abs(value) > tol:
            print(f'{i}, {j}, {k}, {l} -- {value}')


def quartic_tensor_from_operator(op, fundamental_operators,
                                 perm_sym=False):

    # -- Convert fundamental operators back and forth from index

    op_idx_map = get_operator_index_map(fundamental_operators)
    op_idx_set = set(op_idx_map)

    nop = len(fundamental_operators)
    h_quart = np.zeros((nop, nop, nop, nop), dtype=complex)
    
    for term in op:
        op_list, prefactor = term
        if len(op_list) == 4:

            d, t = list(zip(*op_list)) # split in two lists with daggers and tuples resp
            t = [tuple(x) for x in t]

            # check creation/annihilation order
            assert( d == (True, True, False, False) ) 
            
            if all([ x in op_idx_set for x in t ]):
                i, j, k, l = [ op_idx_map.index(x) for x in t ]

                if perm_sym:
                    # -- all pair wise permutations:
                    # -- i <-> j and k <-> l
                    # -- with fermionic permutation signs

                    h_quart[i, j, k, l] = +0.25 * prefactor
                    h_quart[j, i, k, l] = -0.25 * prefactor
                    h_quart[j, i, l, k] = +0.25 * prefactor
                    h_quart[i, j, l, k] = -0.25 * prefactor
 
                else:

                    h_quart[i, j, k, l] = prefactor

    return h_quart


def quartic_permutation_symmetrize(U):

    U = 0.5 * ( U - np.swapaxes(U, 0, 1) )
    U = 0.5 * ( U - np.swapaxes(U, 2, 3) )

    return U


def get_gw_tensor(H_int, fundamental_operators):

    """ Takes a TRIQS operator object and extracts the quartic terms
    and returns the corresponding antisymmetrized quartic tensor in
    vertex index order, i.e., cc+cc+. """
    
    U_abcd = quartic_tensor_from_operator(H_int, fundamental_operators)
    #print('--> tens from op')
    #print_tensor(U_abcd)
    
    # -- Group in c^+cc^+c ( from c^+c^+cc )
    U_abcd = np.ascontiguousarray(np.transpose(U_abcd, (0, 3, 1, 2)))    
    #print('--> transp')
    #print_tensor(U_abcd)

    U_abcd = (U_abcd + np.transpose(U_abcd, (2,3,0,1)))
    #print('--> perm')
    #print_tensor(U_abcd)

    return U_abcd


def hubbard_atom_hf():

    beta = 1.0
    mu = 0.0
    nw = 128
    
    U = 0.1
    t = 1.0
    e_u, e_d = 0.0, 0.0

    H_int = U * n('up',0) * n('do',0)
    H_hop = t * (c_dag('up',0) * c('do',0) + c_dag('do',0) * c('up',0)) + \
        e_u * n('up', 0) + e_d * n('do', 0)
    
    print(f'H_hop =\n{H_hop}')
    print(f'H_int =\n{H_int}')

    U = 1 / np.sqrt(2) * np.array([
        [1,  1],
        [1, -1]])

    np.testing.assert_array_almost_equal(U @ U.T.conj(), np.eye(2)) # Check unitarity
    
    gf_struct = [['up', 1], ['do', 1]]
    fundamental_operators = [c('up', 0), c('do', 0)]
    
    Hp_hop = operator_single_particle_transform(H_hop, U, fundamental_operators)    
    Hp_int = operator_single_particle_transform(H_int, U, fundamental_operators)    
    Hp = Hp_hop + Hp_int

    print('-'*72)
    print(f'Hp_hop = {Hp_hop}')
    print(f'Hp_int = {Hp_int}')
    
    V_aaaa = get_gw_tensor(H_int, fundamental_operators)
    Vp_aaaa = get_gw_tensor(Hp_int, fundamental_operators)

    print('-'*72)
    print('--> V_aaaa')
    print_tensor(V_aaaa)
    print('--> Vp_aaaa')
    print_tensor(Vp_aaaa)

    h_hop  = quadratic_matrix_from_operator(H_hop,  fundamental_operators)
    hp_hop = quadratic_matrix_from_operator(Hp_hop, fundamental_operators)

    print('-'*72)
    print(f'h_hop  =\n{h_hop.real}')
    print(f'hp_hop =\n{hp_hop.real}')
    np.testing.assert_array_almost_equal(U @ h_hop @ U.T.conj(), hp_hop) 
    
    tb_opts = dict(
        units = [(1, 0, 0)],
        orbital_positions = [(0,0,0)] * 2,
        orbital_names = ['up_0', 'do_0'],
        )

    H_r = TBLattice(hopping = {(0,): h_hop,}, **tb_opts)
    Hp_r = TBLattice(hopping = {(0,): hp_hop,}, **tb_opts)

    kmesh = H_r.get_kmesh(n_k=(1, 1, 1))
    
    e_k = H_r.fourier(kmesh)
    ep_k = Hp_r.fourier(kmesh)

    # -- Hartree Fock? solution
    
    hs  = HartreeFockSolver(e_k,  beta, H_int,  gf_struct)
    hsp = HartreeFockSolver(ep_k, beta, Hp_int, gf_struct)

    hs.solve_newton_mu(mu=mu)
    hsp.solve_newton_mu(mu=mu)
    print('-'*72)

    rho_aa  = hs.density_matrix()
    rhop_aa = hsp.density_matrix()

    np.testing.assert_array_almost_equal(U @ rho_aa @ U.T.conj(), rhop_aa)
    
    np.testing.assert_array_almost_equal(
        hs.total_density(), hsp.total_density())

    # -- GW solution

    wmesh = MeshImFreq(beta, 'Fermion', nw)
    g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)    
    gp_wk = lattice_dyson_g0_wk(mu=mu, e_k=ep_k, mesh=wmesh)    

    V_k = Gf(mesh=kmesh, target_shape=[2]*4)
    V_k.data[:] = V_aaaa    

    Vp_k = Gf(mesh=kmesh, target_shape=[2]*4)
    Vp_k.data[:] = Vp_aaaa    

    def get_rho(g_wk):
        wmesh = g_wk.mesh[0]
        g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
        g_w.data[:] = np.sum(g_wk.data, axis=1) / len(kmesh)
        rho = g_w.density().real
        return rho

    rho = get_rho(g_wk)
    rhop = get_rho(gp_wk)

    #print('-'*72)
    #print(f'rho =\n{rho}')
    #print(f'rhop =\n{rhop}')

    np.testing.assert_array_almost_equal(U @ rho @ U.T.conj(), rhop) 
    
    max_iter = 100
    N_old = float('inf')
    
    for gw_iter in range(max_iter):
        
        sigma_h_k = hartree_sigma(V_k, g_wk)
        sigma_f_k = fock_sigma(V_k, g_wk)
        sigma_k = sigma_h_k + sigma_f_k

        sigmap_h_k = hartree_sigma(Vp_k, gp_wk)
        sigmap_f_k = fock_sigma(Vp_k, gp_wk)
        sigmap_k = sigmap_h_k + sigmap_f_k

        s = np.squeeze(sigma_k.data)
        sp = np.squeeze(sigmap_k.data)

        np.testing.assert_array_almost_equal(sp, U @ s @ U.T.conj())
        
        sigma_wk = g_wk.copy()
        sigma_wk.data[:] = sigma_k.data[:, None]

        sigmap_wk = gp_wk.copy()
        sigmap_wk.data[:] = sigmap_k.data[:, None]
        
        g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
        gp_wk = lattice_dyson_g_wk(mu, ep_k, sigmap_wk)

        rho = get_rho(g_wk)
        N = np.sum(np.diag(rho))

        rhop = get_rho(gp_wk)
        Np = np.sum(np.diag(rhop))
        
        #print(N, Np)

        dN = np.abs(N_old - N)
        N_old = N
        
        if dN < 1e-6: break

    print(f'rho =\n{rho}')
    print(f'rhop =\n{rhop}')

    np.testing.assert_array_almost_equal(rho_aa, rho)
    np.testing.assert_array_almost_equal(rhop_aa, rhop)
    
    print(f'sigma =\n{sigma_k.data.real}')
    print(f'sigma_h =\n{sigma_h_k.data.real}')
    print(f'sigma_f =\n{sigma_f_k.data.real}')
    print(f'sigmap =\n{sigmap_k.data.real}')
    print(f'sigmap_h =\n{sigmap_h_k.data.real}')
    print(f'sigmap_f =\n{sigmap_f_k.data.real}')

    np.testing.assert_array_almost_equal(hs.M, np.squeeze(sigma_k.data))
    np.testing.assert_array_almost_equal(hsp.M, np.squeeze(sigmap_k.data))

    np.testing.assert_array_almost_equal(
        np.squeeze(sigmap_f_k.data), np.zeros((2,2)))
    
            
def hubbard_dimer_hf():

    mu = 0.0
    beta = 1.0
    nw = 128

    if True:
        t_u, t_d = 1.0, 0.5
        U, U_01 = 3.0, 1.0
        e_u0, e_d0 = 0.0, 1.0
        e_u1, e_d1 = 0.3, -0.3
    else:
        t_u, t_d = 1.0, 1.0
        U, U_01 = 1.0, 1.0
        e_u0, e_d0 = 0.0, 0.0
        e_u1, e_d1 = 0.0, 0.0
    
    # Dimer with dens-dens interaction and hopping
    # vs. molecular orbitals with transformed interaction
    
    H_int = U * n('up',0) * n('do',0) + U * n('up',1) * n('do',1) + \
        U_01 * (n('up',0) + n('do',0)) * (n('up',1) + n('do',1))
    H_hop = \
        e_u0*n('up',0) + e_d0*n('do',0) + \
        e_u1*n('up',1) + e_d1*n('do',1) + \
        t_u * ( c_dag('up',0) * c('up',1) + c_dag('up',1) * c('up', 0) ) + \
        t_d * ( c_dag('do',0) * c('do',1) + c_dag('do',1) * c('do', 0) )
    H = H_hop + H_int
    
    print('-'*72)
    print(f'H_hop = {H_hop}')
    print(f'H_int = {H_int}')
    #print(f'H = {H}')

    # -- Single particle trasform to molecular orbitals
    
    U = 1 / np.sqrt(2) * np.array([
        [1,  1],
        [1, -1]])

    U = block_diag(U, U)
    np.testing.assert_array_almost_equal(U @ U.T.conj(), np.eye(4)) # Check unitarity
    
    gf_struct = [['up', 2], ['do', 2]]
    fundamental_operators = [c('up', 0), c('up',1), c('do', 0), c('do', 1)]
    
    Hp_hop = operator_single_particle_transform(H_hop, U, fundamental_operators)    
    Hp_int = operator_single_particle_transform(H_int, U, fundamental_operators)    
    Hp = Hp_hop + Hp_int

    print('-'*72)
    print(f'Hp_hop = {Hp_hop}')
    print(f'Hp_int = {Hp_int}')

    V_aaaa  = get_gw_tensor(H_int,  fundamental_operators)
    Vp_aaaa = get_gw_tensor(Hp_int, fundamental_operators)

    print('-'*72)
    print('--> V_aaaa')
    for i,j,k,l in itertools.product(range(4), repeat=4):
        value = V_aaaa[i, j, k, l]
        if np.abs(value) > 1e-9:
            print(f'{i}, {j}, {k}, {l} -- {value}')

    print('-'*72)
    print('--> Vp_aaaa')
    for i,j,k,l in itertools.product(range(4), repeat=4):
        value = Vp_aaaa[i, j, k, l]
        if np.abs(value) > 1e-9:
            print(f'{i}, {j}, {k}, {l} -- {value}')

    #exit()
            
    h_hop  = quadratic_matrix_from_operator(H_hop,  fundamental_operators)
    hp_hop = quadratic_matrix_from_operator(Hp_hop, fundamental_operators)

    print('-'*72)
    print(f'h_hop  =\n{h_hop.real}')
    print(f'hp_hop =\n{hp_hop.real}')
    np.testing.assert_array_almost_equal(U @ h_hop @ U.T.conj(), hp_hop) 
    
    # -- Lattices
    
    tb_opts = dict(
        units = [(1, 0, 0)],
        orbital_positions = [(0,0,0)] * 4,
        orbital_names = ['up_0', 'up_1', 'do_0', 'do_1'],
        )

    H_r = TBLattice(hopping = {(0,): h_hop,}, **tb_opts)
    Hp_r = TBLattice(hopping = {(0,): hp_hop,}, **tb_opts)

    kmesh = H_r.get_kmesh(n_k=(1, 1, 1))
    
    e_k = H_r.fourier(kmesh)
    ep_k = Hp_r.fourier(kmesh)

    # -- Hartree Fock? solution
    
    hs  = HartreeFockSolver(e_k,  beta, H_int,  gf_struct)
    hsp = HartreeFockSolver(ep_k, beta, Hp_int, gf_struct)

    hs.solve_newton_mu(mu=mu)
    hsp.solve_newton_mu(mu=mu)
    print('-'*72)

    rho_aa  = hs.density_matrix()
    rhop_aa = hsp.density_matrix()

    np.testing.assert_array_almost_equal(
        hs.total_density(), hsp.total_density())

    #exit()
    
    np.testing.assert_array_almost_equal(U @ rho_aa @ U.T.conj(), rhop_aa) 

    # -- GW solution

    wmesh = MeshImFreq(beta, 'Fermion', nw)
    g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)    
    gp_wk = lattice_dyson_g0_wk(mu=mu, e_k=ep_k, mesh=wmesh)    

    V_k = Gf(mesh=kmesh, target_shape=[4]*4)
    V_k.data[:] = V_aaaa    

    Vp_k = Gf(mesh=kmesh, target_shape=[4]*4)
    Vp_k.data[:] = Vp_aaaa    

    def get_rho(g_wk):
        wmesh = g_wk.mesh[0]
        g_w = Gf(mesh=wmesh, target_shape=g_wk.target_shape)
        g_w.data[:] = np.sum(g_wk.data, axis=1) / len(kmesh)
        rho = g_w.density().real
        return rho

    rho = get_rho(g_wk)
    rhop = get_rho(gp_wk)

    #print('-'*72)
    #print(f'rho =\n{rho}')
    #print(f'rhop =\n{rhop}')

    np.testing.assert_array_almost_equal(U @ rho @ U.T.conj(), rhop) 
    
    max_iter = 100
    N_old = float('inf')
    
    for gw_iter in range(max_iter):
        
        sigma_h_k = hartree_sigma(V_k, g_wk)
        sigma_f_k = fock_sigma(V_k, g_wk)
        sigma_k = sigma_h_k + sigma_f_k

        sigmap_h_k = hartree_sigma(Vp_k, gp_wk)
        sigmap_f_k = fock_sigma(Vp_k, gp_wk)
        sigmap_k = sigmap_h_k + sigmap_f_k

        s = np.squeeze(sigma_k.data)
        sp = np.squeeze(sigmap_k.data)

        np.testing.assert_array_almost_equal(sp, U @ s @ U.T.conj())
        
        sigma_wk = g_wk.copy()
        sigma_wk.data[:] = sigma_k.data[:, None]

        sigmap_wk = gp_wk.copy()
        sigmap_wk.data[:] = sigmap_k.data[:, None]
        
        g_wk = lattice_dyson_g_wk(mu, e_k, sigma_wk)
        gp_wk = lattice_dyson_g_wk(mu, ep_k, sigmap_wk)

        rho = get_rho(g_wk)
        N = np.sum(np.diag(rho))

        rhop = get_rho(gp_wk)
        Np = np.sum(np.diag(rhop))
        
        #print(N, Np)

        dN = np.abs(N_old - N)
        N_old = N
        
        if dN < 1e-6: break

    print(f'rho =\n{rho}')
    print(f'rhop =\n{rhop}')

    np.testing.assert_array_almost_equal(rho_aa, rho)
    np.testing.assert_array_almost_equal(rhop_aa, rhop)
    
    print(f'sigma =\n{sigma_k.data.real}')
    print(f'sigmap =\n{sigmap_k.data.real}')

    np.testing.assert_array_almost_equal(hs.M, np.squeeze(sigma_k.data))
    np.testing.assert_array_almost_equal(hsp.M, np.squeeze(sigmap_k.data))
    

if __name__ == '__main__':

    hubbard_atom_hf()
    hubbard_dimer_hf()
