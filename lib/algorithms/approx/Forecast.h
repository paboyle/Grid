/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/approx/Forecast.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef INCLUDED_FORECAST_H
#define INCLUDED_FORECAST_H

namespace Grid {

  // Abstract base class.
  // Takes a matrix (Mat), a source (phi), and a vector of Fields (chi)
  // and returns a forecasted solution to the system D*psi = phi (psi).
  template<class Matrix, class Field>
  class Forecast
  {
    public:
      virtual Field operator()(Matrix &Mat, const Field& phi, const std::vector<Field>& chi) = 0;
  };

  // Implementation of Brower et al.'s chronological inverter (arXiv:hep-lat/9509012),
  // used to forecast solutions across poles of the EOFA heatbath.
  //
  // Modified from CPS (cps_pp/src/util/dirac_op/d_op_base/comsrc/minresext.C)
  template<class Matrix, class Field>
  class ChronoForecast : public Forecast<Matrix,Field>
  {
    public:
      Field operator()(Matrix &Mat, const Field& phi, const std::vector<Field>& prev_solns)
      {
        int degree = prev_solns.size();
        Field chi(phi); // forecasted solution

        // Trivial cases
        if(degree == 0){ chi = zero; return chi; }
        else if(degree == 1){ return prev_solns[0]; }

        RealD dot;
        ComplexD xp;
        Field r(phi); // residual
        Field Mv(phi);
        std::vector<Field> v(prev_solns); // orthonormalized previous solutions
        std::vector<Field> MdagMv(degree,phi);

        // Array to hold the matrix elements
        std::vector<std::vector<ComplexD>> G(degree, std::vector<ComplexD>(degree));

        // Solution and source vectors
        std::vector<ComplexD> a(degree);
        std::vector<ComplexD> b(degree);

        // Orthonormalize the vector basis
        for(int i=0; i<degree; i++){
          v[i] *= 1.0/std::sqrt(norm2(v[i]));
          for(int j=i+1; j<degree; j++){ v[j] -= innerProduct(v[i],v[j]) * v[i]; }
        }

        // Perform sparse matrix multiplication and construct rhs
        for(int i=0; i<degree; i++){
          b[i] = innerProduct(v[i],phi);
          Mat.M(v[i],Mv);
          Mat.Mdag(Mv,MdagMv[i]);
          G[i][i] = innerProduct(v[i],MdagMv[i]);
        }

        // Construct the matrix
        for(int j=0; j<degree; j++){
        for(int k=j+1; k<degree; k++){
          G[j][k] = innerProduct(v[j],MdagMv[k]);
          G[k][j] = std::conj(G[j][k]);
        }}

        // Gauss-Jordan elimination with partial pivoting
        for(int i=0; i<degree; i++){

          // Perform partial pivoting
          int k = i;
          for(int j=i+1; j<degree; j++){ if(std::abs(G[j][j]) > std::abs(G[k][k])){ k = j; } }
          if(k != i){
            xp = b[k];
            b[k] = b[i];
            b[i] = xp;
            for(int j=0; j<degree; j++){
              xp = G[k][j];
              G[k][j] = G[i][j];
              G[i][j] = xp;
            }
          }

          // Convert matrix to upper triangular form
          for(int j=i+1; j<degree; j++){
            xp = G[j][i]/G[i][i];
            b[j] -= xp * b[i];
            for(int k=0; k<degree; k++){ G[j][k] -= xp*G[i][k]; }
          }
        }

        // Use Gaussian elimination to solve equations and calculate initial guess
        chi = zero;
        r = phi;
        for(int i=degree-1; i>=0; i--){
          a[i] = 0.0;
          for(int j=i+1; j<degree; j++){ a[i] += G[i][j] * a[j]; }
          a[i] = (b[i]-a[i])/G[i][i];
          chi += a[i]*v[i];
          r -= a[i]*MdagMv[i];
        }

        RealD true_r(0.0);
        ComplexD tmp;
        for(int i=0; i<degree; i++){
          tmp = -b[i];
          for(int j=0; j<degree; j++){ tmp += G[i][j]*a[j]; }
          tmp = std::conj(tmp)*tmp;
          true_r += std::sqrt(tmp.real());
        }

        RealD error = std::sqrt(norm2(r)/norm2(phi));
        std::cout << GridLogMessage << "ChronoForecast: |res|/|src| = " << error << std::endl;

        return chi;
      };
  };

}

#endif
