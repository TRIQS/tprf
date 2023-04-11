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

using namespace nda::clef;

namespace {
placeholder<0> iw;
placeholder<1> w;
placeholder<2> n1;
placeholder<3> n2;
placeholder<4> n;
  
} // namespace

namespace triqs_tprf {

namespace hubbard_atom {
  
  typedef std::complex<double> val_t; 
  typedef gf<imfreq, tensor_valued<4>> temp_1d_t;
  typedef gf<prod<imfreq, imfreq>, tensor_valued<4>> temp_2d_t;

  g_iw_t single_particle_greens_function(int nw, double beta, double U) {
    g_iw_t g_iw{{beta, Fermion, nw}, {1, 1}};
    g_iw[iw] << 1./(iw - U*U / (4 * iw) );
    return g_iw;
  }

  g2_iw_t chi_ph_magnetic(int nw, int nwf, double beta, double U) {

    auto mb = mesh::imfreq{beta, Boson, nw};
    auto mf = mesh::imfreq{beta, Fermion, nwf};

    temp_1d_t C{mb, {1, 1, 1, 1}}, D{mb, {1, 1, 1, 1}};
    temp_2d_t a0{{mb, mf}, {1, 1, 1, 1}}, b0{{mb, mf}, {1, 1, 1, 1}};
    temp_2d_t b1{{mb, mf}, {1, 1, 1, 1}}, b2{{mb, mf}, {1, 1, 1, 1}};
    g2_iw_t chi{{mb, mf, mf}, {1, 1, 1, 1}};

    val_t I(0., 1.);
    
    val_t A = I*U*0.5;
    val_t B = I*U*0.5 * sqrt( abs(3. - exp(beta*U*0.5)) / (1 + exp(beta*U*0.5)) );
    C(w) << 0.;
    C[0] = -beta * U * 0.5 / (1 + exp(-beta * U * 0.5)); // set w=0
    D(w) << U*U*0.25 * (1. + C(w))/(1. - C(w));
    
    val_t B0 = 1.;
    val_t A0 = 1.;
    val_t B1 = 1.;
    val_t B2 = I;

    a0(w, n) << A0 * beta*0.5 * (-n*(n+w) - A*A) /
      ( (-n*n + U*U*0.25) * (-(n+w)*(n+w) + U*U*0.25) );

    b0(w, n) << B0 * beta*0.5 * (-n*(n+w) - B*B) /
      ( (-n*n + U*U*0.25) * (-(n+w)*(n+w) + U*U*0.25) );

    b1(w, n) << B1 * sqrt(U*(1-C(w))) * (-n*(n+w) - D(w)) /
      ( (-n*n + U*U*0.25) * (-(n+w)*(n+w) + U*U*0.25) );

    b2(w, n) << B2 * sqrt(U*U*U*0.25) * sqrt(U*U/(1-C(w)) - w*w) /
      ( (-n*n + U*U*0.25) * (-(n+w)*(n+w) + U*U*0.25) );

    chi(w, n1, n2) << kronecker(n1, n2) * (b0(w, n1) + a0(w, n1))
                    + kronecker(n1, -w - n2) * (b0(w, n1) - a0(w, n1))
                    + b1(w, n1) * b1(w, n2) + b2(w, n1) * b2(w, n2);

    return chi;
  }

  g2_iw_t gamma_ph_magnetic(int nw, int nwf, double beta, double U) {

    auto mb = mesh::imfreq{beta, Boson, nw};
    auto mf = mesh::imfreq{beta, Fermion, nwf};

    temp_2d_t a0{{mb, mf}, {1, 1, 1, 1}}, b0{{mb, mf}, {1, 1, 1, 1}};
    temp_1d_t D{mb, {1, 1, 1, 1}};
    g2_iw_t gamma{{mb, mf, mf}, {1, 1, 1, 1}};

    val_t I(0., 1.);

    int s = 1; // Sign +1 (d,m); -1 (s,t)

    val_t A = I*U*0.5;
    // Modified for cplx sqrt
    val_t B = I*U*0.5 * sqrt( abs(3. - exp(beta*U*0.5)) / (1 + exp(beta*U*0.5)) );
    //std::cout << "3. - exp(beta*U*0.5) = " << 3. - exp(beta*U*0.5) << "\n";
    
    val_t B0 = 1.;
    val_t A0 = 1.;
    val_t B1 = 1.;

    a0(w, n) << beta*0.5*A*A/A0 *
      (-n*n+U*U*0.25)*(-(n+s*w)*(n+s*w)+U*U*0.25) /
      ((-n*(n+s*w)-A*A)*(-n*(n+s*w)));

    b0(w, n) << beta*0.5*B*B/B0 * 
      (-n*n+U*U*0.25)*(-(n+s*w)*(n+s*w)+U*U*0.25) /
      ((-n*(n+s*w)-B*B)*(-n*(n+s*w)));

    for( auto const & w : mb ) {
      val_t sqrtBBww = sqrt(4*B*B - w*w);
      val_t powBBUU1 = std::pow(4.*B*B/(U*U) + 1., 2);
      val_t UU4 = U*U*0.25;

      // NB! The s factor has an additional sign cf Thunstrom expression 
      D[w] = -U*abs(B1)*abs(B1)/(B0*B0) * UU4*(UU4*powBBUU1 - w*w) /
	( U*tan( beta*0.25*(sqrtBBww - I*w) ) / sqrtBBww - s );
    }
    
    /*
    // does not work for the sqrt() calls .. ? FIXME
    D(w) << -U * abs(B1)*abs(B1) / (B0*B0) *
      U*U*0.25 * ( U*U*0.25 * std::pow(4*B*B/(U*U) + 1, 2) - w*w) /
      ( U*tan(beta*0.25 * (sqrt(4*B*B-w*w) - I*w)) / sqrt(4*B*B-w*w) - s );
    */
      
    gamma(w, n1, n2) << kronecker(n1, n2) * (b0(w, n1) + a0(w, n1))
                      + kronecker(n1, -w - n2) * (b0(w, n1) - a0(w, n1))
                      + D(w) / (-n1*(n1+s*w) - B*B) / (-n2*(n2+s*w) - B*B)
                      - B1*B1 / (B0*B0) * U;

    gamma(w, n1, n2) << - 1/beta/beta * gamma(w, n1, n2);
    
    return gamma;
  }
  
} // namespace hubbard_atom

} // namespace triqs_tprf
