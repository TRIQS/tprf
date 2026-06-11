################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
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
from triqs.plot.mpl_interface import plt,oplot
import triqs.gf
from h5 import HDFArchive
import itertools
# ----------------------------------------------------------------------

from triqs.gf import MeshBrillouinZone, Idx, Gf, MeshProduct
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs.operators import n, c, c_dag, Operator

# ----------------------------------------------------------------------

from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections

# ----------------------------------------------------------------------

from triqs.operators.util.op_struct import set_operator_structure

def quadratic_operator_to_matrix(operator, gf_struct):

    from triqs.operators.util.extractors import dict_to_matrix, extract_h_dict

    dct = extract_h_dict(operator)
    matrix = dict_to_matrix(dct, gf_struct)

    norb = matrix.shape[0] // 2

    # Transpose to make spin the fastest running index in the compositie spin orbital index
    matrix = matrix.reshape(2, norb, 2, norb)
    matrix = matrix.transpose(1,0, 3,2)
    matrix = matrix.reshape(2*norb, 2*norb)
    
    return matrix

def trace_channel(tensor,operator1,operator2):
    
    return np.einsum('abcd,ab,cd',tensor,operator1,operator2)

num_orbitals=1
operator_N = +1*sum( c_dag('up',i)*c('up',i)+c_dag('do',i)*c('do',i) for i in range(num_orbitals) )
operator_Nup = +1*sum( c_dag('up',i)*c('up',i) for i in range(num_orbitals) )
operator_Ndo = +1*sum( c_dag('do',i)*c('do',i) for i in range(num_orbitals) )
operator_Sz = +0.5*sum( c_dag('up',i)*c('up',i)-c_dag('do',i)*c('do',i) for i in range(num_orbitals) )
operator_Sx = +0.5*sum( c_dag('up',i)*c('do',i)+c_dag('do',i)*c('up',i) for i in range(num_orbitals) )
operator_Sy = +0.5j*sum( c_dag('up',i)*c('do',i)-c_dag('do',i)*c('up',i) for i in range(num_orbitals) )
orb_names = list(range(num_orbitals))
spin_names = ('up', 'do')
gf_struct = set_operator_structure(spin_names, orb_names, off_diag=True)


font = {'size'   : 8}

import matplotlib
matplotlib.rc('font', **font)

Q = (0,0,0)

nwfs = [30,28,26,24,22,20,18] # Fermionic frequency cut-offs

fig,axs=plt.subplots(nrows= 1, ncols=4, figsize=(20,10))

def load_h5(filename):
    print(f'--> Loading: {filename}')
    with HDFArchive(filename, 'r') as a:
        p = a['p']
    return p


def plot_chi_kw(filename, operator1,operator2,operator3,operator4,setlabel,ii):
    p=load_h5(filename)

    mesh = p.chi_kw_dbse.mesh[1]
    index = 0
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_dbse[Idx(Q),wn], operator1, operator1).real for wn in mesh], 
        marker='x',label=rf'$\chi(Q,\omega)$ DBSE nwf={p.nwf}' if setlabel else '',c='C0', alpha= 1-ii/len(nwfs) )
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_bse[Idx(Q),wn], operator1, operator1).real for wn in mesh], 
        marker='+',label=rf'$\chi(Q,\omega)$ BSE nwf={p.nwf}' if setlabel else '',c='C1', alpha= 1-ii/len(nwfs) )

    # Exact
    def f(wn):
        return (2*p.B)*p.M/( (2*p.B)**2+wn**2)
    axs[index].plot([wn.value.imag for wn in mesh], [ f(wn.value.imag) for wn in mesh], c='black',alpha=0.5, marker='.'   )

    axs[index].set_xlim(-5,5)

    index = 1
    axs[index].set_xlim(-5,5)
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_dbse[Idx(Q),wn],operator1, operator2).real for wn in mesh], 
        marker='x',label=rf'$\chi(Q,\omega)$ DBSE nwf={p.nwf}' if setlabel else '',c='C0', alpha= 1-ii/len(nwfs) )
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_bse[Idx(Q),wn],operator1, operator2).real for wn in mesh], 
        marker='+',label=rf'$\chi(Q,\omega)$ BSE nwf={p.nwf}' if setlabel else '',c='C1', alpha= 1-ii/len(nwfs) )

    # Exact
    def g(wn):
        return p.M*wn/( (2*p.B)**2+wn**2)
    axs[index].plot([wn.value.imag for wn in mesh], [ g(wn.value.imag) for wn in mesh], c='black',alpha=0.5, marker='.'   )



    # Sz
    index = 2
    axs[index].set_xlim(-5,5)
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_dbse[Idx(Q),wn], operator3, operator3).real for wn in mesh], 
        marker='x',label=rf'$\chi(Q,\omega)$ DBSE nwf={p.nwf}' if setlabel else '',c='C0', alpha= 1-ii/len(nwfs) )
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_bse[Idx(Q),wn], operator3, operator3).real for wn in mesh], 
        marker='+',label=rf'$\chi(Q,\omega)$ BSE nwf={p.nwf}' if setlabel else '',c='C1', alpha= 1-ii/len(nwfs) )

    # N
    index = 3
    axs[index].set_xlim(-5,5)
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_dbse[Idx(Q),wn], operator4, operator4).real for wn in mesh], 
        marker='x',label=rf'$\chi(Q,\omega)$ DBSE nwf={p.nwf}' if setlabel else '',c='C0', alpha= 1-ii/len(nwfs) )
    axs[index].plot([wn.value.imag for wn in mesh], [ trace_channel(p.chi_kw_bse[Idx(Q),wn], operator4, operator4).real for wn in mesh], 
        marker='+',label=rf'$\chi(Q,\omega)$ BSE nwf={p.nwf}' if setlabel else '',c='C1', alpha= 1-ii/len(nwfs) )

    return

operator1 = quadratic_operator_to_matrix(operator_Sx, gf_struct)
operator2 = quadratic_operator_to_matrix(operator_Sy, gf_struct)
operator3 = quadratic_operator_to_matrix(operator_Sz, gf_struct)
operator4 = quadratic_operator_to_matrix(operator_N, gf_struct)
setlabel=True
for ii,nwf in enumerate(nwfs):
    print("nwf=",nwf)
    plot_chi_kw(f'data_bse_nwf_{nwf:03d}_nk_008.h5', operator1, operator2, operator3, operator4, setlabel,ii)

axs[-1].legend()
#fig.suptitle(f'component={component}')
fig.tight_layout()
plt.savefig('fig_chiG.svg')
