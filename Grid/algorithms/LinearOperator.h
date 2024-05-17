/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/LinearOperator.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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

/////////////////////////////////////////////////////////////////////////////////////////////
// LinearOperators Take a something and return a something.
/////////////////////////////////////////////////////////////////////////////////////////////
//
// Hopefully linearity is satisfied and the AdjOp is indeed the Hermitian Conjugateugate (transpose if real):
//SBase
//   i)  F(a x + b y) = aF(x) + b F(y).
//  ii)  <x|Op|y> = <y|AdjOp|x>^\ast
//
// Would be fun to have a test linearity & Herm Conj function!
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Field> class LinearOperatorBase {
public:
  // Support for coarsening to a multigrid
  virtual void OpDiag (const Field &in, Field &out) = 0; // Abstract base
  virtual void OpDir  (const Field &in, Field &out,int dir,int disp) = 0; // Abstract base
  virtual void OpDirAll  (const Field &in, std::vector<Field> &out) = 0; // Abstract base

  virtual void Op     (const Field &in, Field &out) = 0; // Abstract base
  virtual void AdjOp  (const Field &in, Field &out) = 0; // Abstract base
  virtual void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2)=0;
  virtual void HermOp(const Field &in, Field &out)=0;
  virtual ~LinearOperatorBase(){};
};


/////////////////////////////////////////////////////////////////////////////////////////////
// By sharing the class for Sparse Matrix across multiple operator wrappers, we can share code
// between RB and non-RB variants. Sparse matrix is like the fermion action def, and then
// the wrappers implement the specialisation of "Op" and "AdjOp" to the cases minimising
// replication of code.
//
// I'm not entirely happy with implementation; to share the Schur code between herm and non-herm
// while still having a "OpAndNorm" in the abstract base I had to implement it in both cases
// with an assert trap in the non-herm. This isn't right; there must be a better C++ way to
// do it, but I fear it required multiple inheritance and mixed in abstract base classes
/////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// Construct herm op from non-herm matrix
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class MdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  MdagMLinearOperator(Matrix &Mat): _Mat(Mat){};

  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    _Mat.Mdiag(in,out);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    _Mat.Mdir(in,out,dir,disp);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    _Mat.MdirAll(in,out);
  };
  void Op     (const Field &in, Field &out){
    _Mat.M(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    _Mat.Mdag(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    _Mat.MdagM(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    _Mat.MdagM(in,out);
  }
};

////////////////////////////////////////////////////////////////////
// Construct herm op and shift it for mgrid smoother
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class ShiftedMdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  RealD _shift;
public:
  ShiftedMdagMLinearOperator(Matrix &Mat,RealD shift): _Mat(Mat), _shift(shift){};
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    _Mat.Mdiag(in,out);
    assert(0);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    _Mat.Mdir(in,out,dir,disp);
    assert(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    assert(0);
  };
  void Op     (const Field &in, Field &out){
    _Mat.M(in,out);
    assert(0);
  }
  void AdjOp     (const Field &in, Field &out){
    _Mat.Mdag(in,out);
    assert(0);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    _Mat.MdagM(in,out);
    out = out + _shift*in;
  }
};

////////////////////////////////////////////////////////////////////
// Create a shifted HermOp
////////////////////////////////////////////////////////////////////
template<class Field>
class ShiftedHermOpLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field> &_Mat;
  RealD _shift;
public:
  ShiftedHermOpLinearOperator(LinearOperatorBase<Field> &Mat,RealD shift): _Mat(Mat), _shift(shift){};
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    assert(0);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    assert(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    assert(0);
  };
  void Op     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    HermOp(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot = innerProduct(in,out);
    n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    _Mat.HermOp(in,out);
    out = out + _shift*in;
  }
};


////////////////////////////////////////////////////////////////////
// Wrap an already herm matrix
////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class HermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  HermitianLinearOperator(Matrix &Mat): _Mat(Mat){};
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    _Mat.Mdiag(in,out);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    _Mat.Mdir(in,out,dir,disp);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    _Mat.MdirAll(in,out);
  };
  void Op     (const Field &in, Field &out){
    _Mat.M(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    _Mat.M(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    HermOp(in,out);
    ComplexD dot= innerProduct(in,out); n1=real(dot);
    n2=norm2(out);
  }
  void HermOp(const Field &in, Field &out){
    _Mat.M(in,out);
  }
};

template<class Matrix,class Field>
class NonHermitianLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
public:
  NonHermitianLinearOperator(Matrix &Mat): _Mat(Mat){};
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    _Mat.Mdiag(in,out);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    _Mat.Mdir(in,out,dir,disp);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    _Mat.MdirAll(in,out);
  };
  void Op     (const Field &in, Field &out){
    _Mat.M(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    _Mat.Mdag(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    assert(0);
  }
  void HermOp(const Field &in, Field &out){
    assert(0);
  }
};

//////////////////////////////////////////////////////////
// Even Odd Schur decomp operators; there are several
// ways to introduce the even odd checkerboarding
//////////////////////////////////////////////////////////

template<class Field>
class SchurOperatorBase :  public LinearOperatorBase<Field> {
 public:
  virtual  void Mpc      (const Field &in, Field &out) =0;
  virtual  void MpcDag   (const Field &in, Field &out) =0;
  virtual  void MpcDagMpc(const Field &in, Field &out) {
    Field tmp(in.Grid());
    tmp.Checkerboard() = in.Checkerboard();
    Mpc(in,tmp);
    MpcDag(tmp,out);
  }
  virtual void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    out.Checkerboard() = in.Checkerboard();
    MpcDagMpc(in,out);
    ComplexD dot= innerProduct(in,out); 
    n1=real(dot);
    n2=norm2(out);
  }
  virtual void HermOp(const Field &in, Field &out){
    out.Checkerboard() = in.Checkerboard();
    MpcDagMpc(in,out);
  }
  void Op     (const Field &in, Field &out){
    Mpc(in,out);
  }
  void AdjOp     (const Field &in, Field &out){ 
    MpcDag(in,out);
  }
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {
    assert(0); // must coarsen the unpreconditioned system
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    assert(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    assert(0);
  };
};
template<class Matrix,class Field>
  class SchurDiagMooeeOperator :  public SchurOperatorBase<Field> {
 public:
    Matrix &_Mat;
    SchurDiagMooeeOperator (Matrix &Mat): _Mat(Mat){};
    virtual  void Mpc      (const Field &in, Field &out) {
      Field tmp(in.Grid());
      tmp.Checkerboard() = !in.Checkerboard();
      
      _Mat.Meooe(in,tmp);
      _Mat.MooeeInv(tmp,out);
      _Mat.Meooe(out,tmp);
      _Mat.Mooee(in,out);
      axpy(out,-1.0,tmp,out);
    }
    virtual void MpcDag   (const Field &in, Field &out){
      Field tmp(in.Grid());
	
      _Mat.MeooeDag(in,tmp);
      _Mat.MooeeInvDag(tmp,out);
      _Mat.MeooeDag(out,tmp);
      _Mat.MooeeDag(in,out);
      axpy(out,-1.0,tmp,out);
    }
};
template<class Matrix,class Field>
  class SchurDiagOneOperator :  public SchurOperatorBase<Field> {
 protected:
    Matrix &_Mat;
 public:
    SchurDiagOneOperator (Matrix &Mat): _Mat(Mat){};
    
    virtual void Mpc      (const Field &in, Field &out) {
      Field tmp(in.Grid());

      _Mat.Meooe(in,out);
      _Mat.MooeeInv(out,tmp);
      _Mat.Meooe(tmp,out);
      _Mat.MooeeInv(out,tmp);
      axpy(out,-1.0,tmp,in);
    }
    virtual void MpcDag   (const Field &in, Field &out){
      Field tmp(in.Grid());
      
      _Mat.MooeeInvDag(in,out);
      _Mat.MeooeDag(out,tmp);
      _Mat.MooeeInvDag(tmp,out);
      _Mat.MeooeDag(out,tmp);
      axpy(out,-1.0,tmp,in);
    }
};
template<class Matrix,class Field>
  class SchurDiagTwoOperator :  public SchurOperatorBase<Field> {
 protected:
    Matrix &_Mat;
 public:
    SchurDiagTwoOperator (Matrix &Mat): _Mat(Mat){};
    
    virtual void Mpc      (const Field &in, Field &out) {
      Field tmp(in.Grid());
      
      _Mat.MooeeInv(in,out);
      _Mat.Meooe(out,tmp);
      _Mat.MooeeInv(tmp,out);
      _Mat.Meooe(out,tmp);
      
      axpy(out,-1.0,tmp,in);
    }
    virtual  void MpcDag   (const Field &in, Field &out){
      Field tmp(in.Grid());

      _Mat.MeooeDag(in,out);
      _Mat.MooeeInvDag(out,tmp);
      _Mat.MeooeDag(tmp,out);
      _Mat.MooeeInvDag(out,tmp);

      axpy(out,-1.0,tmp,in);
    }
};

template<class Field>
class NonHermitianSchurOperatorBase :  public LinearOperatorBase<Field> 
{
 public:
  virtual void  Mpc      (const Field& in, Field& out) = 0;
  virtual void  MpcDag   (const Field& in, Field& out) = 0;
  virtual void  MpcDagMpc(const Field& in, Field& out) {
    Field tmp(in.Grid());
    tmp.Checkerboard() = in.Checkerboard();
    Mpc(in,tmp);
    MpcDag(tmp,out);
  }
  virtual void HermOpAndNorm(const Field& in, Field& out, RealD& n1, RealD& n2) {
    assert(0);
  }
  virtual void HermOp(const Field& in, Field& out) {
    assert(0);
  }
  void Op(const Field& in, Field& out) {
    Mpc(in, out);
  }
  void AdjOp(const Field& in, Field& out) { 
    MpcDag(in, out);
  }
  // Support for coarsening to a multigrid
  void OpDiag(const Field& in, Field& out) {
    assert(0); // must coarsen the unpreconditioned system
  }
  void OpDir(const Field& in, Field& out, int dir, int disp) {
    assert(0);
  }
  void OpDirAll(const Field& in, std::vector<Field>& out){
    assert(0);
  };
};

template<class Matrix, class Field>
class NonHermitianSchurDiagMooeeOperator :  public NonHermitianSchurOperatorBase<Field> 
{
 public:
  Matrix& _Mat;
 NonHermitianSchurDiagMooeeOperator(Matrix& Mat): _Mat(Mat){};
  virtual void Mpc(const Field& in, Field& out) {
    Field tmp(in.Grid());
    tmp.Checkerboard() = !in.Checkerboard();
    
    _Mat.Meooe(in, tmp);
    _Mat.MooeeInv(tmp, out);
    _Mat.Meooe(out, tmp);
    
    _Mat.Mooee(in, out);
    
    axpy(out, -1.0, tmp, out);
  }
  virtual void MpcDag(const Field& in, Field& out) {
    Field tmp(in.Grid());
    
    _Mat.MeooeDag(in, tmp);
    _Mat.MooeeInvDag(tmp, out);
    _Mat.MeooeDag(out, tmp);
	  
    _Mat.MooeeDag(in, out);
    
    axpy(out, -1.0, tmp, out);
  }
};
    
template<class Matrix,class Field>
class NonHermitianSchurDiagOneOperator : public NonHermitianSchurOperatorBase<Field> 
{
 protected:
  Matrix &_Mat;
  
 public:
  NonHermitianSchurDiagOneOperator (Matrix& Mat): _Mat(Mat){};
  virtual void Mpc(const Field& in, Field& out) {
    Field tmp(in.Grid());
	  
    _Mat.Meooe(in, out);
    _Mat.MooeeInv(out, tmp);
    _Mat.Meooe(tmp, out);
    _Mat.MooeeInv(out, tmp);

    axpy(out, -1.0, tmp, in);
  }
  virtual void MpcDag(const Field& in, Field& out) {
    Field tmp(in.Grid());
    
    _Mat.MooeeInvDag(in, out);
    _Mat.MeooeDag(out, tmp);
    _Mat.MooeeInvDag(tmp, out);
    _Mat.MeooeDag(out, tmp);
    
    axpy(out, -1.0, tmp, in);
  }
};

template<class Matrix, class Field>
class NonHermitianSchurDiagTwoOperator : public NonHermitianSchurOperatorBase<Field> 
{
 protected:
  Matrix& _Mat;
  
 public:
 NonHermitianSchurDiagTwoOperator(Matrix& Mat): _Mat(Mat){};

  virtual void Mpc(const Field& in, Field& out) {
    Field tmp(in.Grid());
    
    _Mat.MooeeInv(in, out);
    _Mat.Meooe(out, tmp);
    _Mat.MooeeInv(tmp, out);
    _Mat.Meooe(out, tmp);

    axpy(out, -1.0, tmp, in);
  }
  virtual void MpcDag(const Field& in, Field& out) {
    Field tmp(in.Grid());
    
    _Mat.MeooeDag(in, out);
    _Mat.MooeeInvDag(out, tmp);
    _Mat.MeooeDag(tmp, out);
    _Mat.MooeeInvDag(out, tmp);

    axpy(out, -1.0, tmp, in);
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Left  handed Moo^-1 ; (Moo - Moe Mee^-1 Meo) psi = eta  -->  ( 1 - Moo^-1 Moe Mee^-1 Meo ) psi = Moo^-1 eta
// Right handed Moo^-1 ; (Moo - Moe Mee^-1 Meo) Moo^-1 Moo psi = eta  -->  ( 1 - Moe Mee^-1 Meo Moo^-1) phi=eta ; psi = Moo^-1 phi
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class Matrix,class Field> using SchurDiagOneRH = SchurDiagTwoOperator<Matrix,Field> ;
template<class Matrix,class Field> using SchurDiagOneLH = SchurDiagOneOperator<Matrix,Field> ;
///////////////////////////////////////////////////////////////////////////////////////////////////
//  Staggered use
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class Matrix,class Field>
class SchurStaggeredOperator :  public SchurOperatorBase<Field> {
 protected:
  Matrix &_Mat;
  Field tmp;
  RealD mass;
 public:
  SchurStaggeredOperator (Matrix &Mat): _Mat(Mat), tmp(_Mat.RedBlackGrid()) 
  { 
    assert( _Mat.isTrivialEE() );
    mass = _Mat.Mass();
  }
  virtual void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    Mpc(in,out);
    ComplexD dot= innerProduct(in,out);
    n1 = real(dot);
    n2 =0.0;
  }
  virtual void HermOp(const Field &in, Field &out){
    Mpc(in,out);
    //    _Mat.Meooe(in,out);
    //    _Mat.Meooe(out,tmp);
    //    axpby(out,-1.0,mass*mass,tmp,in);
  }
  virtual  void Mpc      (const Field &in, Field &out) 
  {
    Field tmp(in.Grid());
    Field tmp2(in.Grid());
	
    //    _Mat.Mooee(in,out);
    //    _Mat.Mooee(out,tmp);

    _Mat.Meooe(in,out);
    _Mat.Meooe(out,tmp);
    axpby(out,-1.0,mass*mass,tmp,in);
  }
  virtual  void MpcDag   (const Field &in, Field &out){
    Mpc(in,out);
  }
  virtual void MpcDagMpc(const Field &in, Field &out) {
    assert(0);// Never need with staggered
  }
};
template<class Matrix,class Field> using SchurStagOperator = SchurStaggeredOperator<Matrix,Field>;

/////////////////////////////////////////////////////////////
// Base classes for functions of operators
/////////////////////////////////////////////////////////////
template<class Field> class OperatorFunction {
public:
  virtual void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) = 0;
  virtual void operator() (LinearOperatorBase<Field> &Linop, const std::vector<Field> &in,std::vector<Field> &out) {
    assert(in.size()==out.size());
    for(int k=0;k<in.size();k++){
      (*this)(Linop,in[k],out[k]);
    }
  };
  virtual ~OperatorFunction(){};
};

template<class Field> class LinearFunction {
public:
  virtual void operator() (const Field &in, Field &out) = 0;

  virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
  {
    assert(in.size() == out.size());

    for (unsigned int i = 0; i < in.size(); ++i)
    {
      (*this)(in[i], out[i]);
    }
  }
  virtual ~LinearFunction(){};
};

template<class Field> class IdentityLinearFunction : public LinearFunction<Field> {
public:
  void operator() (const Field &in, Field &out){
    out = in;
  };
};


/////////////////////////////////////////////////////////////
// Base classes for Multishift solvers for operators
/////////////////////////////////////////////////////////////
template<class Field> class OperatorMultiFunction {
public:
  virtual void operator() (LinearOperatorBase<Field> &Linop, const Field &in, std::vector<Field> &out) = 0;
};

// FIXME : To think about

// Chroma functionality list defining LinearOperator
/*
  virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const = 0;
  virtual void operator() (T& chi, const T& psi, enum PlusMinus isign, Real epsilon) const
  virtual const Subset& subset() const = 0;
  virtual unsigned long nFlops() const { return 0; }
  virtual void deriv(P& ds_u, const T& chi, const T& psi, enum PlusMinus isign) const
  class UnprecLinearOperator : public DiffLinearOperator<T,P,Q>
  const Subset& subset() const {return all;}
  };
*/

////////////////////////////////////////////////////////////////////////////////////////////
// Hermitian operator Linear function and operator function
////////////////////////////////////////////////////////////////////////////////////////////
template<class Field>
class HermOpOperatorFunction : public OperatorFunction<Field> {
  void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {
    Linop.HermOp(in,out);
  };
};

template<typename Field>
class PlainHermOp : public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator();
  LinearOperatorBase<Field> &_Linop;
      
  PlainHermOp(LinearOperatorBase<Field>& linop) : _Linop(linop) 
  {}
      
  void operator()(const Field& in, Field& out) {
    _Linop.HermOp(in,out);
  }
};

template<typename Field>
class FunctionHermOp : public LinearFunction<Field> {
public:
  using LinearFunction<Field>::operator(); 
  OperatorFunction<Field>   & _poly;
  LinearOperatorBase<Field> &_Linop;
      
  FunctionHermOp(OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop) 
    : _poly(poly), _Linop(linop) {};
      
  void operator()(const Field& in, Field& out) {
    _poly(_Linop,in,out);
  }
};

template<class Field>
class Polynomial : public OperatorFunction<Field> {
private:
  std::vector<RealD> Coeffs;
public:
  using OperatorFunction<Field>::operator();

  Polynomial(std::vector<RealD> &_Coeffs) : Coeffs(_Coeffs) { };

  // Implement the required interface
  void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) {

    Field AtoN(in.Grid());
    Field Mtmp(in.Grid());
    AtoN = in;
    out = AtoN*Coeffs[0];
    for(int n=1;n<Coeffs.size();n++){
      Mtmp = AtoN;
      Linop.HermOp(Mtmp,AtoN);
      out=out+AtoN*Coeffs[n];
    }
  };
};

NAMESPACE_END(Grid);
