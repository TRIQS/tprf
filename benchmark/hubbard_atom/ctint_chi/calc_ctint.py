# ----------------------------------------------------------------------

from pytriqs.gf import *
from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from triqs_ctint import SolverCore, version

from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import c, c_dag, n
from itertools import product

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def Gf_from_struct(mesh, struct):
    # Without block structure
    if not isinstance(struct[0], list):
        return Gf(mesh=mesh, indices=struct)

    # With block structure
    G_lst = []
    for bl, idx_lst in struct:
        G_lst.append(Gf(mesh=mesh, indices=idx_lst))
    return BlockGf(name_list=[bl[0] for bl in struct], block_list=G_lst)

# ----------------------------------------------------------------------
if __name__ == '__main__':

    p = ParameterCollection(
        beta = 5.,
        U = 5.,
        mu = 2.5,
        h = 0.0,
        n_iw = 100,
        delta = 0.1,
        gf_struct = [ ['up',[0]], ['dn',[0]] ],
        store_list = ['G_iw', 'G2_iw', 'G2c_iw'],
        )

    p.solver_core = ParameterCollection(
        beta = p.beta, 
        gf_struct = dict(p.gf_struct),
        n_iw = p.n_iw,  
        n_tau = 100001
        )

    p.solve = ParameterCollection(
        n_warmup_cycles = int(1e4),
        #n_cycles = 10000000,
        n_cycles = int(1e4) / mpi.size,
        length_cycle = 20,
        measure_M_tau = True,
        measure_M4_iw = True,
        n_iw_M4 = 40,
        n_iW_M4 = 1,
        post_process = True,      
        )

    p.version = ParameterCollection()
    p.version.grab_attribs(version, ['version', 'ctint_hash', 'triqs_hash'])

    # -- The alpha tensor

    p.diag = 0.5 + p.delta
    p.odiag = 0.5 - p.delta
    p.solve.alpha = [
        [[p.diag, p.odiag] for i in indices ] for bl, indices in p.gf_struct ]
    
    # -- Total impurity hamiltonian and interaction part

    p.H_0 = - p.mu*( n('up',0) + n('dn',0) ) - p.h*( n('up',0) - n('dn',0) )
    p.H_int = p.U * n('up',0) * n('dn',0)
    p.H = p.H_0 + p.H_int
    p.solve.h_int = p.H_int

    # -- Non-Interacting Impurity Green function

    iw_mesh = MeshImFreq(p.beta, 'Fermion', p.n_iw)
    p.G0_iw = Gf_from_struct(mesh=iw_mesh, struct=p.gf_struct)
        
    h_field_dict = dict( up=p.h, dn=-p.h )
    for name, g0_iw in p.G0_iw:
        h_field = h_field_dict[name]
        g0_iw << inverse(iOmega_n + p.mu + h_field)
        
    # -- CTINT
    
    S = SolverCore(**p.solver_core.dict())
    S.G0_iw << p.G0_iw
    S.solve(**p.solve.dict())

    # -- Store to hdf5 archive
    
    p.grab_attribs(S, p.store_list)
    
    if mpi.is_master_node():
        with HDFArchive("data_ctint.h5", 'w') as results:
            results["p"] = p
