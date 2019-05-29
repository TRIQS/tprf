
from common import *

from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

def plot_ps(ps):    

    subp = [2, 1, 1]

    plt.subplot(*subp); subp[-1] += 1
    plt.plot(ps.iter, ps.dG_l, 's-', label=r'$\Delta G_l$')
    plt.plot(ps.iter, ps.dM, 'o-', label=r'$\Delta M$')
    plt.semilogy([], [])
    plt.ylabel('$\Delta G_l$, $\Delta M$')
    plt.legend(loc='best')
    plt.xlabel('Iteration')

    plt.subplot(*subp); subp[-1] += 1
    plt.plot(ps.iter, ps.B, 's-', label=r'$B$')
    plt.plot(ps.iter, ps.M, 'o-', label=r'$M$')
    plt.ylabel(r'$M$, $B$')
    plt.legend(loc='best')
    plt.xlabel('Iteration')

    plt.tight_layout()
    plt.savefig('figure_field_sc.svg')
    
def plot_p(p):    

    xlim = [-50, 50]
    
    subp = [3, 1, 1]
    
    plt.subplot(*subp); subp[-1] += 1
    oploti(p.g0_w[0,0], label=r'$G_0(i\omega_n)$')
    oploti(p.g_w[0,0], label=r'$G(i\omega_n)$')
    oploti(p.sigma_w[0,0], label=r'$\Sigma(i\omega_n)$')
    plt.legend(loc='best'); plt.xlim(xlim)

    plt.subplot(*subp); subp[-1] += 1
    oplotr(p.G_tau_raw['up'], alpha=0.75)
    oplotr(p.G_tau['up'])
    plt.gca().legend().set_visible(False)

    plt.subplot(*subp); subp[-1] += 1
    G_l = p.G_l.copy()
    for b, g in G_l: G_l[b].data[:] = np.abs(g.data)
    oplotr(G_l['up'], 'o-')
    
    plt.semilogy([], [])
    plt.gca().legend().set_visible(False)

    plt.tight_layout()
    plt.savefig('figure_field_sc_gf.svg')
    
filenames = np.sort(glob.glob('data_B_*.h5'))
ps = []
for filename in filenames:
    with HDFArchive(filename, 'r') as a:
        ps += a['ps'].objects

ps = ParameterCollections(ps)

with HDFArchive('data_B_0.000000.h5', 'r') as a:
    p = a['ps'].objects[-1]

plt.figure(figsize=(3.25*2, 5))
plot_p(p)

plt.figure(figsize=(3.25*2, 5))
plot_ps(ps)

plt.show()
