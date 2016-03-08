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
#ifndef  GRID_ALGORITHM_LINEAR_OP_H
#define  GRID_ALGORITHM_LINEAR_OP_H

namespace Grid {

  /////////////////////////////////////////////////////////////////////////////////////////////
  // LinearOperators Take a something and return a something.
  /////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Hopefully linearity is satisfied and the AdjOp is indeed the Hermitian conjugateugate (transpose if real):
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

      virtual void Op     (const Field &in, Field &out) = 0; // Abstract base
      virtual void AdjOp  (const Field &in, Field &out) = 0; // Abstract base
      virtual void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2)=0;
      virtual void HermOp(const Field &in, Field &out)=0;
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
      void Op     (const Field &in, Field &out){
	_Mat.M(in,out);
      }
      void AdjOp     (const Field &in, Field &out){
	_Mat.Mdag(in,out);
      }
      void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
	_Mat.MdagM(in,out,n1,n2);
      }
      void HermOp(const Field &in, Field &out){
	RealD n1,n2;
	HermOpAndNorm(in,out,n1,n2);
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
      void Op     (const Field &in, Field &out){
	_Mat.M(in,out);
	assert(0);
      }
      void AdjOp     (const Field &in, Field &out){
	_Mat.Mdag(in,out);
	assert(0);
      }
      void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
	_Mat.MdagM(in,out,n1,n2);
	out = out + _shift*in;

	ComplexD dot;	
	dot= innerProduct(in,out);
	n1=real(dot);
	n2=norm2(out);
      }
      void HermOp(const Field &in, Field &out){
	RealD n1,n2;
	HermOpAndNorm(in,out,n1,n2);
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
      void Op     (const Field &in, Field &out){
	_Mat.M(in,out);
      }
      void AdjOp     (const Field &in, Field &out){
	_Mat.M(in,out);
      }
      void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
	ComplexD dot;

	_Mat.M(in,out);
	
	dot= innerProduct(in,out);
	n1=real(dot);

	dot = innerProduct(out,out);
	n2=real(dot);
      }
      void HermOp(const Field &in, Field &out){
	_Mat.M(in,out);
      }
    };

    //////////////////////////////////////////////////////////
    // Even Odd Schur decomp operators; there are several
    // ways to introduce the even odd checkerboarding
    //////////////////////////////////////////////////////////

    template<class Field>
      class SchurOperatorBase :  public LinearOperatorBase<Field> {
    public:
      virtual  RealD Mpc      (const Field &in, Field &out) =0;
      virtual  RealD MpcDag   (const Field &in, Field &out) =0;
      virtual void MpcDagMpc(const Field &in, Field &out,RealD &ni,RealD &no) {
	Field tmp(in._grid);
	ni=Mpc(in,tmp);
	no=MpcDag(tmp,out);
      }
      void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
	MpcDagMpc(in,out,n1,n2);
      }
      void HermOp(const Field &in, Field &out){
	RealD n1,n2;
	HermOpAndNorm(in,out,n1,n2);
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

    };
    template<class Matrix,class Field>
      class SchurDiagMooeeOperator :  public SchurOperatorBase<Field> {
    protected:
      Matrix &_Mat;
    public:
      SchurDiagMooeeOperator (Matrix &Mat): _Mat(Mat){};
      virtual  RealD Mpc      (const Field &in, Field &out) {
	Field tmp(in._grid);
//	std::cout <<"grid pointers: in._grid="<< in._grid << " out._grid=" << out._grid << "  _Mat.Grid=" << _Mat.Grid() << " _Mat.RedBlackGrid=" << _Mat.RedBlackGrid() << std::endl;

	_Mat.Meooe(in,tmp);
	_Mat.MooeeInv(tmp,out);
	_Mat.Meooe(out,tmp);

	_Mat.Mooee(in,out);
	return axpy_norm(out,-1.0,tmp,out);
      }
      virtual  RealD MpcDag   (const Field &in, Field &out){
	Field tmp(in._grid);

	_Mat.MeooeDag(in,tmp);
	_Mat.MooeeInvDag(tmp,out);
	_Mat.MeooeDag(out,tmp);

	_Mat.MooeeDag(in,out);
	return axpy_norm(out,-1.0,tmp,out);
      }
    };
    template<class Matrix,class Field>
      class SchurDiagOneOperator :  public SchurOperatorBase<Field> {
    protected:
      Matrix &_Mat;
    public:
      SchurDiagOneOperator (Matrix &Mat): _Mat(Mat){};

      virtual  RealD Mpc      (const Field &in, Field &out) {
	Field tmp(in._grid);

	_Mat.Meooe(in,out);
	_Mat.MooeeInv(out,tmp);
	_Mat.Meooe(tmp,out);
	_Mat.MooeeInv(out,tmp);

	return axpy_norm(out,-1.0,tmp,in);
      }
      virtual  RealD MpcDag   (const Field &in, Field &out){
	Field tmp(in._grid);

	_Mat.MooeeInvDag(in,out);
	_Mat.MeooeDag(out,tmp);
	_Mat.MooeeInvDag(tmp,out);
	_Mat.MeooeDag(out,tmp);

	return axpy_norm(out,-1.0,tmp,in);
      }
    };

    template<class Matrix,class Field>
      class SchurDiagTwoOperator :  public SchurOperatorBase<Field> {
    protected:
      Matrix &_Mat;
    public:
      SchurDiagTwoOperator (Matrix &Mat): _Mat(Mat){};

      virtual  RealD Mpc      (const Field &in, Field &out) {
	Field tmp(in._grid);

	_Mat.MooeeInv(in,out);
	_Mat.Meooe(out,tmp);
	_Mat.MooeeInv(tmp,out);
	_Mat.Meooe(out,tmp);

	return axpy_norm(out,-1.0,tmp,in);
      }
      virtual  RealD MpcDag   (const Field &in, Field &out){
	Field tmp(in._grid);

	_Mat.MeooeDag(in,out);
	_Mat.MooeeInvDag(out,tmp);
	_Mat.MeooeDag(tmp,out);
	_Mat.MooeeInvDag(out,tmp);

	return axpy_norm(out,-1.0,tmp,in);
      }
    };


    /////////////////////////////////////////////////////////////
    // Base classes for functions of operators
    /////////////////////////////////////////////////////////////
    template<class Field> class OperatorFunction {
    public:
      virtual void operator() (LinearOperatorBase<Field> &Linop, const Field &in, Field &out) = 0;
    };

    template<class Field> class LinearFunction {
    public:
      virtual void operator() (const Field &in, Field &out) = 0;
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


}

#endif
