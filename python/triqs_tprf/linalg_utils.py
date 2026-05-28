# -*- coding: utf-8 -*-

################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2024, Hugo U. R. Strand
# Authors: H. U.R. Strand
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

from math import isqrt

from triqs.gfs import Gf


class ChannelOrder():
    
    """ This class (re-)defines the channel pairing of single-particle indices
    in the three two-particle response channels:
    
    1. particle-hole (PH),
    2. particle-hole-bar (PH_bar), and
    3. particle-particle (PP).

    (For the generalization including fermionic frequencies see channel_grouping.hpp)
    
    The single particle labels are indexed according to:
    
    {a, b, c, d} <=> {0, 1, 2, 3}

    Channel::PH
    -----------

    in the particle-hole channel (PH) the indices are grouped as
    
    {a, b}, {d, c} <=> {0, 1}, {3, 2}

    Channel::PH_bar
    ---------------

    in the particle-hole-bar channel (PH_bar) the indices are grouped as
    
    {a, d}, {c, b} <=> {0, 3}, {2, 1}

    Channel::PP
    -----------

    in the particle-particle channel (Channel_t::PP) the indices are grouped as
    
    {a, c}, {b, d} <=> {0, 2}, {1, 3}
    
    """

    PH     = np.array([0, 1, 3, 2])
    PH_bar = np.array([0, 3, 2, 1])
    PP     = np.array([0, 2, 1, 3])


def matrix_from_tensor(tensor, channel_order=ChannelOrder.PH):

    ts = np.array(tensor.shape)
    n = ts[0]
    
    assert( len(ts) == 4 )
    assert( (n == ts).all() )
    
    return np.transpose(tensor, channel_order).reshape([n**2]*2)


def tensor_from_matrix(matrix, channel_order=ChannelOrder.PH):

    ms = np.array(matrix.shape)
    n2 = ms[0]
    n = isqrt(n2)

    assert( len(ms) == 2 )
    assert( (n2 == ms).all() )
    assert( n**2 == n2 )
    
    return np.transpose(matrix.reshape([n]*4), channel_order)

    
def gf_matrix_from_tensor(G_tensor, channel_order=ChannelOrder.PH):
    
    ts = np.array(G_tensor.target_shape)
    n = ts[0]

    assert( len(ts) == 4 )
    assert( (n == ts).all() )
        
    tensor_data_shape = G_tensor.data.shape

    mesh_shape = tensor_data_shape[:-4] # Strip last 4 tensor dimensions
    matrix_data_shape = np.concatenate((mesh_shape, [n**2]*2))

    nm = len(mesh_shape) # Number of meshes

    # Permute last 4 indices according to channel grouping
    channel_permuted_axes = np.concatenate((np.arange(nm), nm + channel_order))
        
    G_matrix = Gf(mesh=G_tensor.mesh, target_shape=[n**2]*2)

    # Transpose and reshape
    G_matrix.data[:] = np.transpose(G_tensor.data, axes=channel_permuted_axes).reshape(matrix_data_shape)

    return G_matrix


def gf_tensor_from_matrix(G_matrix, channel_order=ChannelOrder.PH):

    ms = np.array(G_matrix.target_shape)
    n2 = ms[0]
    n = isqrt(n2)

    assert( len(ms) == 2 )
    assert( (n2 == ms).all() )
    assert( n**2 == n2 )

    matrix_data_shape = G_matrix.data.shape

    mesh_shape = matrix_data_shape[:-2] # Strip last 2 matrix dimensions
    tensor_data_shape = np.concatenate((mesh_shape, [n]*4))

    nm = len(mesh_shape) # Number of meshes

    # Permute last 4 indices according to channel grouping
    channel_permuted_axes = np.concatenate((np.arange(nm), nm + channel_order))

    G_tensor = Gf(mesh=G_matrix.mesh, target_shape=[n]*4)
    
    # Transpose and reshape
    G_tensor.data[:] = np.transpose(G_matrix.data.reshape(tensor_data_shape), axes=channel_permuted_axes)

    return G_tensor
