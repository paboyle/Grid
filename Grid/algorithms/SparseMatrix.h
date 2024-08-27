/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/SparseMatrix.h

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
#ifndef  GRID_ALGORITHM_SPARSE_MATRIX_H
#define  GRID_ALGORITHM_SPARSE_MATRIX_H


NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////////////////////////////////////////////////////
// Interface defining what I expect of a general sparse matrix, such as a Fermion action
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Field> class SparseMatrixBase {
public:
  virtual GridBase *Grid(void) =0;
  // Full checkerboar operations
  virtual void  M    (const Field &in, Field &out)=0;
  virtual void  Mdag (const Field &in, Field &out)=0;
  virtual void  MdagM(const Field &in, Field &out) {
    Field tmp (in.Grid());
    M(in,tmp);
    Mdag(tmp,out);
  }
  virtual void  MMdag(const Field &in, Field &out) {
    Field tmp (in.Grid());
    Mdag(in,tmp);
    M(tmp,out);
  }
  virtual  void Mdiag    (const Field &in, Field &out)=0;
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp)=0;
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)=0;
  virtual ~SparseMatrixBase() {};
};

/////////////////////////////////////////////////////////////////////////////////////////////
// Interface augmented by a red black sparse matrix, such as a Fermion action
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Field> class CheckerBoardedSparseMatrixBase : public SparseMatrixBase<Field> {
public:
  virtual GridBase *RedBlackGrid(void)=0;

  //////////////////////////////////////////////////////////////////////
  // Query the even even properties to make algorithmic decisions
  //////////////////////////////////////////////////////////////////////
  virtual RealD  Mass(void)        { return 0.0; };
  virtual int    ConstEE(void)     { return 1; }; // Disable assumptions unless overridden
  virtual int    isTrivialEE(void) { return 0; }; // by a derived class that knows better

  // half checkerboard operaions
  virtual  void Meooe    (const Field &in, Field &out)=0;
  virtual  void Mooee    (const Field &in, Field &out)=0;
  virtual  void MooeeInv (const Field &in, Field &out)=0;

  virtual  void MeooeDag    (const Field &in, Field &out)=0;
  virtual  void MooeeDag    (const Field &in, Field &out)=0;
  virtual  void MooeeInvDag (const Field &in, Field &out)=0;
  virtual ~CheckerBoardedSparseMatrixBase() {};
};

NAMESPACE_END(Grid);

#endif
