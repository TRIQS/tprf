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

from triqs.gfs import Gf
from triqs.mesh import MeshImFreq, MeshProduct


from triqs_tprf.linalg_utils import ChannelOrder
from triqs_tprf.linalg_utils import matrix_from_tensor, tensor_from_matrix
from triqs_tprf.linalg_utils import gf_matrix_from_tensor, gf_tensor_from_matrix


def test_conversion_between_tensor_and_matrix_valued_channel_representations():

    channel_orders = [ChannelOrder.PH, ChannelOrder.PH_bar, ChannelOrder.PP]
    
    beta = 1.0
    mf = MeshImFreq(beta, 'Fermion', n_iw=3)
    mb = MeshImFreq(beta, 'Boson', n_iw=10)
    n = 4

    for channel_order in channel_orders:
        for i in range(10):
            U = np.random.random((n, n, n, n))
            U_mat = matrix_from_tensor(U, channel_order=channel_order)
            U_ref = tensor_from_matrix(U_mat, channel_order=channel_order)
            np.testing.assert_array_almost_equal(U, U_ref)
    
    Pi_w = Gf(mesh=mb, target_shape=[n]*4)
    Pi_wnn = Gf(mesh=MeshProduct(mb, mf, mf), target_shape=[n]*4)

    for channel_order in channel_orders:
        for i in range(10):

            Pi_w.data[:] = np.random.random(Pi_w.data.shape)
            Pi_w_mat = gf_matrix_from_tensor(Pi_w, channel_order=channel_order)
            Pi_w_ref = gf_tensor_from_matrix(Pi_w_mat, channel_order=channel_order)

            np.testing.assert_array_almost_equal(Pi_w.data, Pi_w_ref.data)

            Pi_wnn.data[:] = np.random.random(Pi_wnn.data.shape)
            Pi_wnn_mat = gf_matrix_from_tensor(Pi_wnn, channel_order=channel_order)
            Pi_wnn_ref = gf_tensor_from_matrix(Pi_wnn_mat, channel_order=channel_order)

            np.testing.assert_array_almost_equal(Pi_wnn.data, Pi_wnn_ref.data)

    
    # RPA example using TRIQS internal inverse
    
    # This only works for Gf with single MeshImFreq mesh
    # - Constants are treated like unit matrix times constant
    # - Products are matrix products

    from triqs.gfs import inverse

    U = np.random.random((n, n, n, n)) # Need channel grouping for this as well
    U_mat = matrix_from_tensor(U)

    Pi_w_mat = gf_matrix_from_tensor(Pi_w, channel_order=ChannelOrder.PH)
    
    Chi_w_mat = Pi_w_mat.copy()
    Chi_w_mat << Pi_w_mat * inverse(1. - U_mat * Pi_w_mat)

    Chi_w = gf_tensor_from_matrix(Chi_w_mat, channel_order=ChannelOrder.PH)

    # Triangle vertex RPA equation
    
    L_wn = Gf(mesh=MeshProduct(mb, mf), target_shape=[n]*4)
    L_wn_mat = gf_matrix_from_tensor(L_wn, channel_order=ChannelOrder.PH)

    Lirr_wn_mat = L_wn_mat.copy()

    for inu in mf:
        Lirr_wn_mat[:, inu] << L_wn_mat[:, inu] * inverse(1 - U_mat * L_wn_mat[:, inu])
        


if __name__ == "__main__":
    test_conversion_between_tensor_and_matrix_valued_channel_representations()
