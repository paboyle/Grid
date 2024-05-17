/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/NormalEquations.h

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
#ifndef GRID_NORMAL_EQUATIONS_H
#define GRID_NORMAL_EQUATIONS_H

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Take a matrix and form an NE solver calling a Herm solver
///////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Field> class NormalEquations : public LinearFunction<Field>{
private:
  SparseMatrixBase<Field> & _Matrix;
  OperatorFunction<Field> & _HermitianSolver;
  LinearFunction<Field>   & _Guess;
public:

  /////////////////////////////////////////////////////
  // Wrap the usual normal equations trick
  /////////////////////////////////////////////////////
 NormalEquations(SparseMatrixBase<Field> &Matrix, OperatorFunction<Field> &HermitianSolver,
		 LinearFunction<Field> &Guess) 
   :  _Matrix(Matrix), _HermitianSolver(HermitianSolver), _Guess(Guess) {}; 

  void operator() (const Field &in, Field &out){
 
    Field src(in.Grid());
    Field tmp(in.Grid());

    MdagMLinearOperator<SparseMatrixBase<Field>,Field> MdagMOp(_Matrix);
    _Matrix.Mdag(in,src);
    _Guess(src,out);
    _HermitianSolver(MdagMOp,src,out);  // Mdag M out = Mdag in

  }     
};

template<class Field> class HPDSolver : public LinearFunction<Field> {
private:
  LinearOperatorBase<Field> & _Matrix;
  OperatorFunction<Field> & _HermitianSolver;
  LinearFunction<Field>   & _Guess;
public:

  /////////////////////////////////////////////////////
  // Wrap the usual normal equations trick
  /////////////////////////////////////////////////////
 HPDSolver(LinearOperatorBase<Field> &Matrix,
	   OperatorFunction<Field> &HermitianSolver,
	   LinearFunction<Field> &Guess) 
   :  _Matrix(Matrix), _HermitianSolver(HermitianSolver), _Guess(Guess) {}; 

  void operator() (const Field &in, Field &out){
 
    _Guess(in,out);
    _HermitianSolver(_Matrix,in,out);  //M out = in

  }     
};


template<class Field> class MdagMSolver : public LinearFunction<Field> {
private:
  SparseMatrixBase<Field> & _Matrix;
  OperatorFunction<Field> & _HermitianSolver;
  LinearFunction<Field>   & _Guess;
public:

  /////////////////////////////////////////////////////
  // Wrap the usual normal equations trick
  /////////////////////////////////////////////////////
 MdagMSolver(SparseMatrixBase<Field> &Matrix, OperatorFunction<Field> &HermitianSolver,
	     LinearFunction<Field> &Guess) 
   :  _Matrix(Matrix), _HermitianSolver(HermitianSolver), _Guess(Guess) {}; 

  void operator() (const Field &in, Field &out){
 
    MdagMLinearOperator<SparseMatrixBase<Field>,Field> MdagMOp(_Matrix);
    _Guess(in,out);

    _HermitianSolver(MdagMOp,in,out);  // Mdag M out = Mdag in

  }     
};

NAMESPACE_END(Grid);
#endif
