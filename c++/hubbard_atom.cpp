/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "hubbard_atom.hpp"

#include <triqs/clef.hpp>
using namespace triqs::clef;

namespace {
placeholder<0> iw;
placeholder<1> w;
placeholder<2> n1;
placeholder<3> n2;

placeholder<4> a;
placeholder<5> b;
placeholder<6> c;
placeholder<7> d;
} // namespace

namespace tprf {

namespace hubbard_atom {
  
  g_iw_t single_particle_greens_function(int nw, double beta, double U) {
    g_iw_t g_iw{{beta, Fermion, nw}, {1, 1}};
    g_iw[iw] << 1./(iw - U*U / (4 * iw) );
    return g_iw;
  }

  g2_iw_t chi_ph_magnetic(int nw, int nwf, double beta, double U) {

    typedef std::complex<double> val_t; 
    typedef gf<imfreq, tensor_valued<4>> temp_1d_t;
    typedef gf<cartesian_product<imfreq, imfreq>, tensor_valued<4>> temp_2d_t;

    auto mb = gf_mesh<imfreq>{beta, Boson, nw};
    auto mf = gf_mesh<imfreq>{beta, Fermion, nwf};

    temp_1d_t C{mb, {1, 1, 1, 1}};
    temp_1d_t D{mb, {1, 1, 1, 1}};
    
    temp_2d_t a0{{mb, mf}, {1, 1, 1, 1}};
    temp_2d_t b0{{mb, mf}, {1, 1, 1, 1}};
    temp_2d_t b1{{mb, mf}, {1, 1, 1, 1}};
    temp_2d_t b2{{mb, mf}, {1, 1, 1, 1}};

    g2_iw_t chi{{mb, mf, mf}, {1, 1, 1, 1}};
    
    C(w) << 0.;
    C[{mb, 0}] = -beta*U*0.5 / ( 1 + exp(-beta*U*0.5) ); // set zeroth Matsubara frequency    
    D(w) << U*U*0.25 * (1. + C(w))/(1. + C(w));

    val_t A0 = 1.;
    val_t A = 1.j * U * 0.5;

    val_t B0 = 1.;
    val_t B = U*0.5 * sqrt( (3. - exp(beta*U*0.5)) / (1 + exp(beta*U*0.5)) );	
    
    val_t B1 = 1.;
    val_t B2 = 1.j;
    
    for (auto const &w : mb) {
      for (auto const &n : mf) {
	a0[n, w] = A0 * beta*0.5 * (n*(n+w) - A*A) / ((n*n + U*U*0.25) * ((n+w)*(n+w) + U*U * 0.25) );
	b0[n, w] = B0 * beta*0.5 * (n*(n+w) - B*B) / ((n*n + U*U*0.25) * ((n+w)*(n+w) + U*U * 0.25) );
	b1[n, w] = B1 * sqrt(U*(1-C[w])) * (n*(n+w) - D[w]) / ( (n*n + U*U*0.25) * ((n+w)*(n+w) + U*U * 0.25) );
	b2[n, w] = B2 * sqrt(U*U*U*0.25) * sqrt(U*U/(1 - C[w]) + w*w) / ( (n*n + U*U*0.25) * ((n+w)*(n+w) + U*U * 0.25) );
      }
    }

    for (auto const &w : mb) {
      for (auto const &n1 : mf) {
	for (auto const &n2 : mf) {
	  if(n1.index() == n2.index()) chi[w, n1, n2] = b0[w, n1] + a0[w, n1];
	  if(n1.index() == (- w - n2).n) chi[w, n1, n2] = b0[w, n1] - a0[w, n1];
	  chi[w, n1, n2] = b1[n1, w] * b1[n2, w] + b2[n1, w] * b2[n2, w];
	}
      }
    }

    /*
    chi(w, n1, n2) << kronecker(n1, n2) * (b0(w, n1) + a0(w, n1));
                    + kronecker(n1, -w - n2) * (b0(w, n1) - a0(w, n1))
                    + b1[n1, w] * b1[n2, w] + b2[n1, w] * b2[n2, w];
    */
    
    return chi;
  }
  
} // namespace hubbard_atom

} // namespace tprf
