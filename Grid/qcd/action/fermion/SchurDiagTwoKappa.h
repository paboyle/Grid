/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: SchurDiagTwoKappa.h

    Copyright (C) 2017

Author: Christoph Lehner
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#pragma once

NAMESPACE_BEGIN(Grid);

// This is specific to (Z)mobius fermions
template<class Matrix, class Field>
class KappaSimilarityTransform {
public:
  INHERIT_IMPL_TYPES(Matrix);
  Vector<Coeff_t> kappa, kappaDag, kappaInv, kappaInvDag;

  KappaSimilarityTransform (Matrix &zmob) {
    for (int i=0;i<(int)zmob.bs.size();i++) {
      Coeff_t k = 1.0 / ( 2.0 * (zmob.bs[i] *(4 - zmob.M5) + 1.0) );
      kappa.push_back( k );
      kappaDag.push_back( conj(k) );
      kappaInv.push_back( 1.0 / k );
      kappaInvDag.push_back( 1.0 / conj(k) );
    }
  }

  template<typename vobj>
  void sscale(const Lattice<vobj>& in, Lattice<vobj>& out, Coeff_t* s) {
    GridBase *grid=out.Grid();
    out.Checkerboard() = in.Checkerboard();
    assert(grid->_simd_layout[0] == 1); // should be fine for ZMobius for now
    int Ls = grid->_rdimensions[0];
    thread_for(ss, grid->oSites(),
    {
      vobj tmp = s[ss % Ls]*in[ss];
      vstream(out[ss],tmp);
    });
  }

  RealD sscale_norm(const Field& in, Field& out, Coeff_t* s) {
    sscale(in,out,s);
    return norm2(out);
  }

  virtual RealD M       (const Field& in, Field& out) { return sscale_norm(in,out,&kappa[0]);   }
  virtual RealD MDag    (const Field& in, Field& out) { return sscale_norm(in,out,&kappaDag[0]);}
  virtual RealD MInv    (const Field& in, Field& out) { return sscale_norm(in,out,&kappaInv[0]);}
  virtual RealD MInvDag (const Field& in, Field& out) { return sscale_norm(in,out,&kappaInvDag[0]);}

};

template<class Matrix,class Field>
class SchurDiagTwoKappaOperator :  public SchurOperatorBase<Field> {
public:
  KappaSimilarityTransform<Matrix, Field> _S;
  SchurDiagTwoOperator<Matrix, Field> _Mat;

  SchurDiagTwoKappaOperator (Matrix &Mat): _S(Mat), _Mat(Mat) {};

  virtual  RealD Mpc      (const Field &in, Field &out) {
    Field tmp(in.Grid());

    _S.MInv(in,out);
    _Mat.Mpc(out,tmp);
    return _S.M(tmp,out);

  }
  virtual  RealD MpcDag   (const Field &in, Field &out){
    Field tmp(in.Grid());

    _S.MDag(in,out);
    _Mat.MpcDag(out,tmp);
    return _S.MInvDag(tmp,out);
  }
};

NAMESPACE_END(Grid);


