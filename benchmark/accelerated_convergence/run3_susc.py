from common import *

from triqs.operators import n

if mpi.is_master_node():
    with HDFArchive('data_B_0.000000.h5', 'r') as a:
        p = a['ps'].objects[-1]
else: p = None
p = mpi.bcast(p)

# -- Sample chi_sz

Sz = n('up',0)-n('do',0)

p.solve.n_cycles = int(1e9)
p.solve.measure_G2_iw_ph = False
p.solve.measure_O_tau=(Sz,Sz)
measure_O_tau_min_ins=100

cthyb = triqs_cthyb.Solver(**p.init.dict())
cthyb.G0_iw << p.G0_w
cthyb.solve(**p.solve.dict())
p.chi_tau = cthyb.O_tau

if mpi.is_master_node():
    with HDFArchive('data_chi.h5', 'w') as a: a['p'] = p
