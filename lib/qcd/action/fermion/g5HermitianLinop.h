    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/g5HermitianLinop.h

    Copyright (C) 2015

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
#ifndef G5_HERMITIAN_LINOP
#define G5_HERMITIAN_LINOP

namespace Grid {
  namespace QCD {

////////////////////////////////////////////////////////////////////
// Wrap an already herm matrix
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class Gamma5R5HermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  Gamma5R5HermitianLinearOperator(Matrix &Mat): _Mat(Mat){};
  void Op     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void OpDiag (const Field &in, Field &out) {
    Field tmp(in._grid);
    _Mat.Mdiag(in,tmp);
    G5R5(out,tmp);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    Field tmp(in._grid);
    _Mat.Mdir(in,tmp,dir,disp);
    G5R5(out,tmp);
  }

  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){

    HermOp(in,out);
    
    ComplexD dot;
    dot= innerProduct(in,out);
    n1=real(dot);
    
    dot = innerProduct(out,out);
    n2=real(dot);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in._grid);
    _Mat.M(in,tmp);
    G5R5(out,tmp);
  }
};


template<class Matrix,class Field>
class Gamma5HermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Gamma g5;
public:
    Gamma5HermitianLinearOperator(Matrix &Mat): _Mat(Mat), g5(Gamma::Gamma5) {};
  void Op     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void OpDiag (const Field &in, Field &out) {
    Field tmp(in._grid);
    _Mat.Mdiag(in,tmp);
    out=g5*tmp;
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    Field tmp(in._grid);
    _Mat.Mdir(in,tmp,dir,disp);
    out=g5*tmp;
  }

  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){

    HermOp(in,out);
    
    ComplexD dot;
    dot= innerProduct(in,out);
    n1=real(dot);
    
    dot = innerProduct(out,out);
    n2=real(dot);
  }
  void HermOp(const Field &in, Field &out){
    Field tmp(in._grid);
    _Mat.M(in,tmp);
    out=g5*tmp;
  }
};


}}
#endif
