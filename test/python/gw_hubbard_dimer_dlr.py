################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2023 by Hugo U.R. Strand
# Author: H. U.R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


import numpy as np


from triqs.gf import Gf, MeshImFreq, Idx
from triqs.gf import inverse, iOmega_n
from triqs.lattice.tight_binding import TBLattice


from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import lattice_dyson_g_wk

from triqs_tprf.lattice import fourier_wk_to_wr
from triqs_tprf.lattice import chi_wr_from_chi_wk

from triqs_tprf.gw_solver import GWSolver


class GWHubbardDimer:

    def __init__(
            self,
            beta=20.0, U=1.5, t=1.0, mu=0.0, nw=1024, maxiter=100,
            self_interaction=False, spinless=False,
            gw_flag=True, hartree_flag=False, fock_flag=False):
        
        wmesh = MeshImFreq(beta, 'Fermion', nw)

        if spinless:
            tb_opts = dict(
                units = [(1, 0, 0)],
                orbital_positions = [(0,0,0)],
                orbital_names = ['0'],
                )
            I = np.eye(1)
            
        else:
            tb_opts = dict(
                units = [(1, 0, 0)],
                orbital_positions = [(0,0,0)] * 2,
                orbital_names = ['up_0', 'do_0'],
                )
            I = np.eye(2)

        # Have to use t/2 hopping in the TBLattice since it considers
        # hoppings in both directions, which doubles the total hopping
        # for the Hubbard dimer.

        H_r = TBLattice(hopping = {
            (+1,): -0.5 * t * I,
            (-1,): -0.5 * t * I,
            }, **tb_opts)

        kmesh = H_r.get_kmesh(n_k=(2, 1, 1))
        self.e_k = H_r.fourier(kmesh)

        if self_interaction:
            
            V_aaaa = np.zeros((2, 2, 2, 2))

            V_aaaa[0, 0, 0, 0] = U
            V_aaaa[1, 1, 1, 1] = U
            
            V_aaaa[1, 1, 0, 0] = U
            V_aaaa[0, 0, 1, 1] = U

            self.V_aaaa = V_aaaa
            
        if spinless:

            V_aaaa = np.zeros((1, 1, 1, 1))
            V_aaaa[0, 0, 0, 0] = U
            self.V_aaaa = V_aaaa
            
        if not spinless and not self_interaction:
            
            from triqs.operators import n, c, c_dag, Operator, dagger
            from triqs_tprf.gw import get_gw_tensor
            
            self.H_int = U * n('up',0) * n('do',0)
            self.fundamental_operators = [c('up', 0), c('do', 0)]
            self.V_aaaa = get_gw_tensor(self.H_int, self.fundamental_operators)
        
        self.V_k = Gf(mesh=kmesh, target_shape=self.V_aaaa.shape)
        self.V_k.data[:] = self.V_aaaa    

        gw = GWSolver(self.e_k, self.V_k, wmesh, mu=mu)
        gw.solve_iter(
            maxiter=maxiter,
            gw=gw_flag, hartree=hartree_flag, fock=fock_flag,
            spinless=spinless)
        gw.calc_real_space()
        
        self.gw = gw

        for key, val in gw.__dict__.items():
            setattr(self, key, val)


def test_gw_hubbard_dimer(verbose=False):
 
    """
    Comparing to analytical expressions from:
    Chapter 4: Hubbard Dimer in GW and Beyond, by Pina Romaniello

    In the book:
    Simulating Correlations with Computers - Modeling and Simulation Vol. 11
    E. Pavarini and E. Koch (eds.)
    Forschungszentrum Ju ̈lich, 2021, ISBN 978-3-95806-529-1
    
    https://www.cond-mat.de/events/correl21/manuscripts/correl21.pdf    
    """
    
    beta = 20.0
    U = 1.5
    t = 1.0
    nw = 1024
    mu = 0.0
    
    gw = GWHubbardDimer(
        beta = beta,
        U = U,
        t = t,
        mu = mu,
        nw = 1024,
        maxiter = 1,
        spinless = True,
        hartree_flag = True,
        )

    gw_tensor = GWHubbardDimer(
        beta = beta,
        U = U,
        t = t,
        mu = mu,
        nw = 1024,
        maxiter = 1,
        self_interaction = True, # Use tensor structure with V_0000 = U
        hartree_flag = False,
        )

    np.testing.assert_array_almost_equal(
        gw.P_wk[0,0,0,0].data, gw_tensor.P_wk[0,0,0,0].data)

    np.testing.assert_array_almost_equal(
        gw.W_wk[0,0,0,0].data, gw_tensor.W_wk[0,0,0,0].data)

    # The self-interaction calc, does not produce correct Hartree sigma, fix here..
    gw_tensor.sigma_wk.data[:] += U/2
    
    np.testing.assert_array_almost_equal(
        gw.sigma_wk[0,0].data, gw_tensor.sigma_wk[0,0].data)
    
    wmesh = gw.g_wk.mesh[0]
    bmesh = gw.P_wk.mesh[0]
    kmesh = gw.P_wk.mesh[1]

    W_wk_ref = gw.W_wk.copy()

    for k in kmesh:
        for w in bmesh:

            # -- The extra factor of 2 in the denominator
            # -- is given by the interaction with two types of fermions.
            # (one of those contributions is strictly speaking self-interaction)

            # Analytic result Eq. (30)
            
            W_wk_ref[w, k][0,0,0,0] = U/(1 - 2*U * gw.P_wk[w, k][0,0,0,0])

    W_wr_ref = chi_wr_from_chi_wk(W_wk_ref)

    np.testing.assert_array_almost_equal(W_wr_ref.data, gw.W_wr.data)

    # -- Print index structure of the relevant quantities
    
    print(f'g_wr0 =\n{gw.g_wr[Idx(0), Idx(0, 0, 0)]}')
    print(f'g_wr1 =\n{gw.g_wr[Idx(0), Idx(1, 0, 0)]}')

    print('V_aaaa =')
    print_tensor(gw.V_aaaa)

    print('V_wr0')
    print_tensor(gw.V_wr[Idx(0), Idx(0, 0, 0)])
    print('V_wr1')
    print_tensor(gw.V_wr[Idx(0), Idx(1, 0, 0)])

    print('P_wr0')
    print_tensor(gw.P_wr[Idx(0), Idx(0, 0, 0)])
    print('P_wr1')
    print_tensor(gw.P_wr[Idx(0), Idx(1, 0, 0)])

    print('W_wr0')
    print_tensor(gw.W_wr[Idx(0), Idx(0, 0, 0)])
    print('W_wr1')
    print_tensor(gw.W_wr[Idx(0), Idx(1, 0, 0)])

    i,j,k,l = 0,0,0,0

    # == Compare GW result with analytic expression
    # == from the book chapter in the doc string above.
    
    # -- Analytic expression Eq. (18) for G_0

    g0_0_w = Gf(mesh=wmesh, target_shape=[])
    g0_1_w = Gf(mesh=wmesh, target_shape=[])
    
    for w in wmesh:
        g0_0_w[w] = +0.5/(w - t) + 0.5/(w + t)
        g0_1_w[w] = -0.5/(w - t) + 0.5/(w + t)

    np.testing.assert_array_almost_equal(
        gw.g0_wr[:, Idx(0, 0, 0)][0, 0].data, g0_0_w.data)

    np.testing.assert_array_almost_equal(
        gw.g0_wr[:, Idx(1, 0, 0)][0, 0].data, g0_1_w.data)
    
    # -- Analytic expression Eq. (29) for P

    P_0_w = Gf(mesh=bmesh, target_shape=[])
    P_1_w = Gf(mesh=bmesh, target_shape=[])

    for w in bmesh:
        P_0_w[w] = + 0.25 / (w - 2*t) - 0.25 / (w + 2*t)
        P_1_w[w] = - 0.25 / (w - 2*t) + 0.25 / (w + 2*t)
        
    np.testing.assert_array_almost_equal(
        gw.P_wr[:, Idx(0, 0, 0)][0,0,0,0].data, P_0_w.data)

    np.testing.assert_array_almost_equal(
        gw.P_wr[:, Idx(1, 0, 0)][0,0,0,0].data, P_1_w.data)

    # -- Analytic expression Eq. (30) for W

    W_0_w = Gf(mesh=bmesh, target_shape=[])
    W_1_w = Gf(mesh=bmesh, target_shape=[])

    e = 0.0
    h2 = 4*t**2 + 4*U*t # Eq. in text after Eq. (31)    
    h = np.sqrt(h2)

    for w in bmesh:
        W_0_w[w] = U + 2 * U**2 * t / (complex(w)**2 - h2)
        W_1_w[w] = 0 - 2 * U**2 * t / (complex(w)**2 - h2)
        
    np.testing.assert_array_almost_equal(
        gw.W_wr[:, Idx(0, 0, 0)][i,j,k,l].data, W_0_w.data)

    np.testing.assert_array_almost_equal(
        gw.W_wr[:, Idx(1, 0, 0)][i,j,k,l].data, W_1_w.data)

    # Analytic expression in Eq. (31) for Sigma = G_0 W_0
    
    sigma_0_w = Gf(mesh=wmesh, target_shape=[])
    sigma_1_w = Gf(mesh=wmesh, target_shape=[])

    for w in wmesh:
        sigma_0_w[w] = U/2 + U**2 * t / (2*h) * \
            ( 1/(w - (e + t + h)) + 1/(w - (e - t - h)) )
        sigma_1_w[w] = U**2 * t / (2*h) * \
            ( 1/(w - (e + t + h)) - 1/(w - (e - t - h)) )
        
    np.testing.assert_array_almost_equal(
        gw.sigma_wr[:, Idx(0, 0, 0)][0, 0].data, sigma_0_w.data)

    np.testing.assert_array_almost_equal(
        gw.sigma_wr[:, Idx(1, 0, 0)][0, 0].data, sigma_1_w.data)
        
    # -- Analytic expression Eq. (32) for G = 1/[1/G_0 - Sigma] with Sigma = G_0 W_0
    
    g_0_w = Gf(mesh=wmesh, target_shape=[])
    g_1_w = Gf(mesh=wmesh, target_shape=[])

    A = np.sqrt((2*t + h + U/2)**2 + 4*U**2*t/h)
    B = np.sqrt((2*t + h - U/2)**2 + 4*U**2*t/h)

    w1_p = 0.5 * (2*e - h + U/2 + A)
    w1_m = 0.5 * (2*e - h + U/2 - A)

    w2_p = 0.5 * (2*e + h + U/2 + B)
    w2_m = 0.5 * (2*e + h + U/2 - B)

    R1 = (h + 2*t + U/2) / (4 * A)
    R2 = (-h -2*t + U/2) / (4 * B)
    
    for w in wmesh:
        g_0_w[w] = \
            + (0.25 + R1)/(w - w1_p) + (0.25 - R1)/(w - w1_m) + \
            + (0.25 + R2)/(w - w2_p) + (0.25 - R2)/(w - w2_m)
        g_1_w[w] = \
            - (0.25 + R1)/(w - w1_p) - (0.25 - R1)/(w - w1_m) + \
            + (0.25 + R2)/(w - w2_p) + (0.25 - R2)/(w - w2_m)

    np.testing.assert_array_almost_equal(
        gw.g_wr[:, Idx(0, 0, 0)][0, 0].data, g_0_w.data)

    np.testing.assert_array_almost_equal(
        gw.g_wr[:, Idx(1, 0, 0)][0, 0].data, g_1_w.data)

    # -- Visualize the result
    
    if verbose:
        from triqs.plot.mpl_interface import oplot, oploti, oplotr, plt

        plt.figure(figsize=(9, 10))

        subp = [5, 2, 1]
        xlim = [-10, 10]

        def plot_cf(g_tprf, g_ref, subp, plot=oploti):
            plt.subplot(*subp); subp[-1] += 1
            plot(g_ref, 'g.', label='ref')
            plot(g_tprf, 'c-', label='tprf')
            plt.xlim(xlim)


        plot_cf(gw.g0_wr[:, Idx(0, 0, 0)][0, 0], g0_0_w, subp, plot=oploti)
        plt.ylabel(r'$G_0(r=0)$')

        plot_cf(gw.g0_wr[:, Idx(1, 0, 0)][0, 0], g0_1_w, subp, plot=oplotr)
        plt.ylabel(r'$G_0(r=1)$')

        plot_cf(gw.P_wr[:, Idx(0, 0, 0)][0,0,0,0], P_0_w, subp, plot=oplotr)
        plt.ylabel(r'$P(r=0)$')

        plot_cf(gw.P_wr[:, Idx(1, 0, 0)][0,0,0,0], P_1_w, subp, plot=oplotr)
        plt.ylabel(r'$P(r=1)$')

        i,j,k,l = 0,0,0,0

        plot_cf(gw.W_wr[:, Idx(0, 0, 0)][i,j,k,l], W_0_w, subp, plot=oplot)
        plt.ylabel(r'$W(r=0)$')            
        plt.title(f'{i},{j},{k},{l}')

        plot_cf(gw.W_wr[:, Idx(1, 0, 0)][i,j,k,l], W_1_w, subp, plot=oplot)
        plt.ylabel(r'$W(r=1)$')
        plt.title(f'{i},{j},{k},{l}')

        plot_cf(gw.sigma_wr[:, Idx(0, 0, 0)][0, 0], sigma_0_w, subp, plot=oplot)
        plt.ylabel(r'$\Sigma(r=0)$')

        plot_cf(gw.sigma_wr[:, Idx(1, 0, 0)][0, 0], sigma_1_w, subp, plot=oplot)
        plt.ylabel(r'$\Sigma(r=1)$')
            

        plt.subplot(*subp); subp[-1] += 1
        oplot(g_0_w, 'g.', label='ref')
        oplot(gw.g_wr[:, Idx(0, 0, 0)][0, 0], 'c-', label='tprf g')
        oploti(gw.g0_wr[:, Idx(0, 0, 0)][0, 0], 'r-', label='tprf g0')
        plt.xlim(xlim)
        plt.ylabel(r'$G(r=0)$')
            
        plt.subplot(*subp); subp[-1] += 1
        oplot(g_1_w, 'g.')
        oplot(gw.g_wr[:, Idx(1, 0, 0)][0, 0], 'c-', label='tprf g')
        oplotr(gw.g0_wr[:, Idx(1, 0, 0)][0, 0], 'r-', label='tprf g0')
        plt.xlim(xlim)
        plt.ylabel(r'$G(r=1)$')

        
        plt.tight_layout()
        plt.show()


def print_tensor(U, tol=1e-9):
    assert( len(U.shape) == 4)
    n = U.shape[0]

    import itertools
    for i,j,k,l in itertools.product(range(n), repeat=4):
        value = U[i, j, k, l]
        if np.abs(value) > tol:
            print(f'{i}, {j}, {k}, {l} -- {value}')


if __name__ == '__main__':

    test_gw_hubbard_dimer(verbose=True)
