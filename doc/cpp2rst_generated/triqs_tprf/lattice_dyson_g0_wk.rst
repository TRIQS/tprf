..
   Generated automatically by cpp2rst

.. highlight:: c


.. _lattice_dyson_g0_wk:

lattice_dyson_g0_wk
===================

**Synopsis**:

.. code-block:: c

    triqs_tprf::g_wk_t lattice_dyson_g0_wk (double mu, triqs_tprf::e_k_cvt e_k,
   gf_mesh<triqs::gfs::imfreq> mesh)

Construct a non-interacting Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`


Parameters
----------

 * **mu**: chemical potential :math:`\mu`

 * **e_k**: discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`

 * **mesh**: imaginary frequency mesh



Return value
------------

Matsubara frequency lattice Green's function $G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})$

Documentation
-------------


  Computes

  .. math::
     G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
         (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
	 \right]^{-1}_{a\bar{b}},

  using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, chemical potential :math:`\mu`,
  and a Matsubara frequency Green's function mesh.