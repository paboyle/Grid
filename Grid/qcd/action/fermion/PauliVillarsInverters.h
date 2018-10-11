    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/SchurRedBlack.h

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
#pragma once

namespace Grid {
namespace QCD {

template<class Field>
class PauliVillarsSolverUnprec
{
 public:
  ConjugateGradient<Field> & CG;
  PauliVillarsSolverUnprec(  ConjugateGradient<Field> &_CG) : CG(_CG){};

  template<class Matrix>
  void operator() (Matrix &_Matrix,const Field &src,Field &sol)
  {
    RealD m = _Matrix.Mass();
    Field A  (_Matrix.FermionGrid());

    MdagMLinearOperator<Matrix,Field> HermOp(_Matrix);

    _Matrix.SetMass(1.0);
    _Matrix.Mdag(src,A);
    CG(HermOp,A,sol);
    _Matrix.SetMass(m);
  };
};

template<class Field,class SchurSolverType>
class PauliVillarsSolverRBprec
{
 public:
  SchurSolverType & SchurSolver;
  PauliVillarsSolverRBprec( SchurSolverType &_SchurSolver) : SchurSolver(_SchurSolver){};

  template<class Matrix>
  void operator() (Matrix &_Matrix,const Field &src,Field &sol)
  {
    RealD m = _Matrix.Mass();
    Field A  (_Matrix.FermionGrid());

    _Matrix.SetMass(1.0);
    SchurSolver(_Matrix,src,sol);
    _Matrix.SetMass(m);
  };
};

template<class Field,class GaugeField>
class PauliVillarsSolverFourierAccel
{
 public:
  GaugeField      & Umu;
  ConjugateGradient<Field> & CG;

  PauliVillarsSolverFourierAccel(GaugeField &_Umu,ConjugateGradient<Field> &_CG) :  Umu(_Umu), CG(_CG)
  {
  };

  template<class Matrix>
  void operator() (Matrix &_Matrix,const Field &src,Field &sol)
  {
    FourierAcceleratedPV<Field, Matrix, typename Matrix::GaugeField > faPV(_Matrix,Umu,CG) ;
    faPV.pvInv(src,sol);
  };
};


}
}
