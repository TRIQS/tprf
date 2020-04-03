# ----------------------------------------------------------------------

import os
import glob
import numpy as np

# ----------------------------------------------------------------------
from pytriqs.gf import *
from h5 import HDFArchive

# ----------------------------------------------------------------------    

from pytriqs.operators import n
from pytriqs.operators.util.op_struct import set_operator_structure

from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt     

# ----------------------------------------------------------------------

from brew_dmft.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
if __name__ == '__main__':

    U = 5.0

    plt.figure(figsize=(3.25, 2*3))
    plt.title('Static Magnetic Susceptibility' + '\n' +
              r'half-filled Hubbard atom $U=%2.2f$' % U)
    
    # ------------------------------------------------------------------
    # -- Analytic result

    T = np.logspace(-3, 2, num=400)
    beta = 1./T

    # -- Analytic magnetization expecation value
    # -- and static susceptibility
    
    Z = 2. + 2*np.exp(-beta * 0.5 * U)
    m2 = 0.25 * (2 / Z)
    chi_m = 2. * beta * m2
    
    plt.plot(1./beta, chi_m, '-k', lw=0.5)

    # ------------------------------------------------------------------
    # -- external field pyed
    
    filenames = glob.glob('pyed_beta*/data_pyed_extrap*.h5')

    style = 'sk'
    for filename in filenames:
        print '--> Loading:', filename

        with HDFArchive(filename, 'r') as s:
            field = s['field']

        plt.plot(1./field.beta, field.chi, style, alpha=0.25)
        plt.plot(1./field.beta, field.chi_exp, '.r', alpha=0.25)

    plt.plot([], [], style, alpha=0.25,
             label=r'ED $\partial_h \langle \hat{m} \rangle$')
    plt.plot([], [], '.r', alpha=0.25,
             label=r'ED $2\beta \langle \hat{m}^2 \rangle$')

    # ------------------------------------------------------------------
    # -- dynamic pomerol

    filenames = glob.glob('dynamic_*/*.h5')

    styles = { 256:'^b', 128:'>r', 64:'<g', 32:'.m', 16:'.c', 8:'.y' }

    for filename in filenames:
        print '--> Loading:', filename
        with HDFArchive(filename, 'r') as s:
            p = s['p']

        if hasattr(p, 'chi'):
            style = styles[p.nwf]
            plt.plot(1./p.beta, p.chi, style, alpha=0.5)

    for nw in np.sort(styles.keys()):
        style = styles[nw]
        plt.plot([], [], style, alpha=0.5,
                 label=r'Analytic Tr$[\chi_m]/\beta^2$, $n_w=%i$' % nw)

    # ------------------------------------------------------------------

    plt.legend(loc='best', fontsize=6, frameon=False)
    #plt.xlabel(r'$\beta$')
    plt.xlabel(r'$T$')
    plt.ylabel(r'$\chi_m$')
    plt.loglog([], [])

    plt.tight_layout()
    plt.savefig('figure_hubbard_atom_static_magnetic_susceptibility.pdf')
    plt.show()
    
