# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation for a Hubbard atom.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com
s
"""

# ----------------------------------------------------------------------

from pytriqs.operators import c, c_dag
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from pomerol2triqs import PomerolED

# ----------------------------------------------------------------------
def make_calc():
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    params = dict(
        beta = 2.0,
        U = 5.0,
        mu = 2.5, # only working for half filling mu = U/2
        ntau = 3,
        niw = 4,
        niw_gf = 4,
        nw_small=2,
        nwf_small=3,
        )

    # ------------------------------------------------------------------

    class Dummy():
        def __init__(self):
            pass

    d = Dummy() # storage space
    d.params = params

    print '--> Solving SIAM with parameters'
    for key, value in params.items():
        print '%10s = %-10s' % (key, str(value))
        globals()[key] = value # populate global namespace
    
    # ------------------------------------------------------------------

    up, do = 'up', 'dn'
    docc = c_dag(up,0) * c(up,0) * c_dag(do,0) * c(do,0)
    nA = c_dag(up,0) * c(up,0) + c_dag(do,0) * c(do,0)
    d.H = -mu * nA + U * docc
    
    # ------------------------------------------------------------------
    # -- Exact diagonalization

    # Conversion from TRIQS to Pomerol notation for operator indices
    # TRIQS:   block_name, inner_index
    # Pomerol: site_label, orbital_index, spin_name
    index_converter = {
        (up, 0) : ('loc', 0, 'up'),
        (do, 0) : ('loc', 0, 'down'),
        }

    # -- Create Exact Diagonalization instance
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(d.H) # -- Diagonalize H

    gf_struct = {up : [0], do : [0]}

    # -- Single-particle Green's functions
    d.G_iw = ed.G_iw(gf_struct, beta, n_iw=niw_gf)
    #d.G_tau = ed.G_tau(gf_struct, beta, n_tau=ntau)
    #d.G_w = ed.G_w(gf_struct, beta, energy_window=(-2.5, 2.5), n_w=100, im_shift=0.01)

    # -- Particle-particle two-particle Matsubara frequency Green's function
    opt = dict(
        beta=beta, gf_struct=gf_struct,
        #blocks=set([("up", "up"), ("up", "dn"), ("dn", "up"), ("dn", "dn")]),
        blocks=set([("up", "up"), ("up", "dn")]),
        n_iw=niw, n_inu=niw)
    
    d.G2_iw_AABB = ed.G2_iw_inu_inup(
        channel='AllFermionic', block_order='AABB', **opt)

    #d.G2_iw_ABBA = ed.G2_iw_inu_inup(
    #    channel='AllFermionic', block_order='ABBA', **opt)

    opt['n_iw'] = nw_small
    opt['n_inu'] = nwf_small
    
    #d.G2_iw_pp = ed.G2_iw_inu_inup(channel='PP', **opt)
    d.G2_iw_ph = ed.G2_iw_inu_inup(channel='PH', **opt)
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_pomerol_hubbard_atom.h5'
    with HDFArchive(filename,'w') as res:
        for key, value in d.__dict__.items():
            res[key] = value
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc()
