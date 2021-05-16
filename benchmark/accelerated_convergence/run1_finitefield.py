from common import *

p = ParameterCollection(
    t=1., B=0., U=10., mu=0., n_k=16, n_iter=10, G_l_tol=2e-5,
    solve = ParameterCollection(
        length_cycle = 10, n_warmup_cycles = 1000, n_cycles = int(2e6),
        move_double = False, measure_G_l = True
        ),
    init = ParameterCollection(
        beta = 1, n_l = 10, n_iw = 400, n_tau = 4000,
        gf_struct = [('up',[0]), ('do',[0])]
        ),
    )

p = setup_dmft_calculation(p)
for B in [0., 0.05, 0.1, 0.15,]:
    p.B = B
    ps = solve_self_consistent_dmft(p)
    p = ps[-1]
    if mpi.is_master_node():
        with HDFArchive('data_B_{:f}.h5'.format(p.B), 'w') as a:
            a['ps'] = ParameterCollections(ps)
