/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacian.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef COVARIANT_LAPLACIAN_H
#define COVARIANT_LAPLACIAN_H

namespace Grid {
namespace QCD {

template <class Impl>
class LaplacianAdjointField {
 public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef SU<Nc>::LatticeAlgebraVector AVector;

  LaplacianAdjointField(GridBase* grid) : U(Nd, grid){};

  void ImportGauge(const GaugeField& _U) {
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
    }
  }

  void Mdiag(const GaugeLinkField& in, GaugeLinkField& out) { assert(0); }

  void Mdir(const GaugeLinkField& in, GaugeLinkField& out, int dir, int disp) { assert(0); }

/*
  // Operator with algebra vector inputs and outputs
  void M2(const AVector& in, AVector& out) {
    double kappa = 0.9;
    //Reconstruct matrix

    GaugeLinkField tmp(in._grid);
    GaugeLinkField tmp2(in._grid);
    GaugeLinkField sum(in._grid);
    GaugeLinkField out_mat(in._grid);
    GaugeLinkField in_mat(in._grid);
    SU<Nc>::FundamentalLieAlgebraMatrix(in, in_mat);

    sum = zero;
    for (int mu = 0; mu < Nd; mu++) {
      tmp  = U[mu] * Cshift(in_mat, mu, +1) * adj(U[mu]);
      tmp2 = adj(U[mu]) * in_mat * U[mu];
      sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_mat;
    }
    out_mat = (1.0 - kappa) * in_mat - kappa/(double(4*Nd)) * sum;
    // Project
    SU<Nc>::projectOnAlgebra(out, out_mat);
  }
*/

    void M(const GaugeLinkField& in, GaugeLinkField& out) {
    double kappa = 0.999;
    //Reconstruct matrix

    GaugeLinkField tmp(in._grid);
    GaugeLinkField tmp2(in._grid);
    GaugeLinkField sum(in._grid);
    sum = zero;
    for (int mu = 0; mu < Nd; mu++) {
      tmp  = U[mu] * Cshift(in, mu, +1) * adj(U[mu]);
      tmp2 = adj(U[mu]) * in * U[mu];
      sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in;
    }
    out = (1.0 - kappa) * in - kappa/(double(4*Nd)) * sum;
  }

 private:
  std::vector<GaugeLinkField> U;
};


template <class Impl>
class LaplacianAlgebraField {
 public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef SU<Nc>::LatticeAlgebraVector AVector;

  LaplacianAlgebraField(GridBase* grid) : U(Nd, grid){};

  void ImportGauge(const GaugeField& _U) {
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
    }
  }

  void Mdiag(const AVector& in, AVector& out) { assert(0); }

  void Mdir(const AVector& in, AVector& out, int dir, int disp) { assert(0); }

  // Operator with algebra vector inputs and outputs
  void M(const AVector& in, AVector& out) {
    double kappa = 0.999;
    //Reconstruct matrix

    GaugeLinkField tmp(in._grid);
    GaugeLinkField tmp2(in._grid);
    GaugeLinkField sum(in._grid);
    GaugeLinkField out_mat(in._grid);
    GaugeLinkField in_mat(in._grid);
    SU<Nc>::FundamentalLieAlgebraMatrix(in, in_mat);

    sum = zero;
    for (int mu = 0; mu < Nd; mu++) {
      tmp  = U[mu] * Cshift(in_mat, mu, +1) * adj(U[mu]);
      tmp2 = adj(U[mu]) * in_mat * U[mu];
      sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_mat;
    }
    out_mat = (1.0 - kappa) * in_mat - kappa/(double(4*Nd)) * sum;
    // Project
    SU<Nc>::projectOnAlgebra(out, out_mat);
  }

 private:
  std::vector<GaugeLinkField> U;
};


}
}

#endif
