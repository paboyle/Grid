/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: 

Copyright (C) 2015-2016

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#pragma once

NAMESPACE_BEGIN(Grid); 

template<class FermionOperatorD, class FermionOperatorF, class SchurOperatorD, class  SchurOperatorF> 
class MixedPrecisionConjugateGradientOperatorFunction : public OperatorFunction<typename FermionOperatorD::FermionField> {

 public:
    typedef typename FermionOperatorD::FermionField FieldD;
    typedef typename FermionOperatorF::FermionField FieldF;

    using OperatorFunction<FieldD>::operator();

    RealD   Tolerance;
    RealD   InnerTolerance; //Initial tolerance for inner CG. Defaults to Tolerance but can be changed
    Integer MaxInnerIterations;
    Integer MaxOuterIterations;
    GridBase* SinglePrecGrid;
    RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance

    FermionOperatorF &FermOpF;
    FermionOperatorD &FermOpD;;
    SchurOperatorF &LinOpF;
    SchurOperatorD &LinOpD;

    Integer TotalInnerIterations; //Number of inner CG iterations
    Integer TotalOuterIterations; //Number of restarts
    Integer TotalFinalStepIterations; //Number of CG iterations in final patch-up step

     MixedPrecisionConjugateGradientOperatorFunction(RealD tol, RealD tolInner,
						    Integer maxinnerit, 
						    Integer maxouterit,
						    GridBase *_SinglePrecGrid,
                                                    FermionOperatorF &_FermOpF,
                                                    FermionOperatorD &_FermOpD,
						    SchurOperatorF   &_LinOpF,
						    SchurOperatorD   &_LinOpD) : 
      LinOpF(_LinOpF),
      LinOpD(_LinOpD),
      FermOpF(_FermOpF),
      FermOpD(_FermOpD),
      Tolerance(tol), 
      InnerTolerance(tolInner), 
      MaxInnerIterations(maxinnerit), 
      MaxOuterIterations(maxouterit), 
      SinglePrecGrid(_SinglePrecGrid),
      OuterLoopNormMult(100.) 
  { assert(tolInner<0.01);    };

    void operator()(LinearOperatorBase<FieldD> &LinOpU, const FieldD &src, FieldD &psi)
    {

      SchurOperatorD * SchurOpU = static_cast<SchurOperatorD *>(&LinOpU);
      
      // Assumption made in code to extract gauge field
      // We could avoid storing LinopD reference alltogether ?
      assert(&(SchurOpU->_Mat)==&(LinOpD._Mat));

      ////////////////////////////////////////////////////////////////////////////////////
      // Moving this to a Clone method of fermion operator would allow to duplicate the 
      // physics parameters and decrease gauge field copies
      ////////////////////////////////////////////////////////////////////////////////////
      auto &Umu_d = FermOpD.GetDoubledGaugeField();
      auto &Umu_f = FermOpF.GetDoubledGaugeField();
      auto &Umu_fe= FermOpF.GetDoubledGaugeFieldE();
      auto &Umu_fo= FermOpF.GetDoubledGaugeFieldO();
      precisionChange(Umu_f,Umu_d);
      pickCheckerboard(Even,Umu_fe,Umu_f);
      pickCheckerboard(Odd ,Umu_fo,Umu_f);

      //////////////////////////////////////////////////////////////////////////////////////////
      // Make a mixed precision conjugate gradient
      //////////////////////////////////////////////////////////////////////////////////////////
      // Could assume red black solver here and remove the SinglePrecGrid parameter???
      MixedPrecisionConjugateGradient<FieldD,FieldF> MPCG(Tolerance, InnerTolerance,MaxInnerIterations,MaxOuterIterations,SinglePrecGrid,LinOpF,LinOpD);

      std::cout << GridLogMessage << "Calling mixed precision Conjugate Gradient src "<<norm2(src) <<std::endl;
      psi=Zero();
      MPCG(src,psi);
    }
  };
NAMESPACE_END(Grid); 
