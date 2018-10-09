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


template<class Field> class Reconstruct5DfromPhysical {
 private:
  OperatorFunction<Field> & _PauliVillarsSolver;
 public:

 /////////////////////////////////////////////////////
 // Wrap the usual normal equations Schur trick
 /////////////////////////////////////////////////////
 Reconstruct5DfromPhysical(OperatorFunction<Field> &PauliVillarsSolver)  :
  _PauliVillarsSolver(PauliVillarsSolver) 
  { };
  /*
  void SliceDump(const Field &f)
  {
    std::vector<TComplex>  C1;
    Field ff(f._grid);
    Gamma G5 ( Gamma::Algebra::Gamma5);
    ff= f+ G5*f;
    ff=ff*0.5;
    {
      auto ip = localInnerProduct(ff,ff);
      sliceSum(ip,C1,0);
      for(int s=0;s<C1.size();s++){
	std::cout << " P+ C[" <<s<<"] = "<<C1[s]<<std::endl;
      }
    }

    ff= f- G5*f;
    ff=ff*0.5;
    {
      auto ip = localInnerProduct(ff,ff);
      sliceSum(ip,C1,0);
      for(int s=0;s<C1.size();s++){
	std::cout << " P- C[" <<s<<"] = "<<C1[s]<<std::endl;
      }
    }

  }
  */

  template<class Matrix>
  void operator() (Matrix & _Matrix,const Field &sol4,const Field &src4, Field &sol5){

    int Ls =  _Matrix.Ls;
    RealD m = _Matrix.Mass();

    Field psi4(_Matrix.GaugeGrid());
    Field psi(_Matrix.FermionGrid());
    Field A  (_Matrix.FermionGrid());
    Field B  (_Matrix.FermionGrid());
    Field c  (_Matrix.FermionGrid());
    Field b  (_Matrix.FermionGrid());

    typedef typename Matrix::Coeff_t Coeff_t;
    MdagMLinearOperator<Matrix,Field> HermOp(_Matrix);

    ///////////////////////////////////////
    //Import source, include Dminus factors
    ///////////////////////////////////////
    _Matrix.ImportPhysicalFermionSource(src4,b); // Includes D_- factor

    ///////////////////////////////////////
    // Set up c from src4
    ///////////////////////////////////////
    std::cout << GridLogMessage<< " ************************************************" << std::endl;
    std::cout << GridLogMessage<< " Reconstrucing 5D propagator using Pauli Villars" << std::endl;
    std::cout << GridLogMessage<< " ************************************************" << std::endl;

    _Matrix.SetMass(1.0); // PauliVillars mass

    _Matrix.Mdag(b,B);   _PauliVillarsSolver(HermOp,B,A);
    _Matrix.Pdag(A,c);      

    _Matrix.SetMass(m);   // Back to physical mass

    //////////////////////////////////////
    // Build Pdag PV^-1 Dm P [-sol4,c2,c3... cL]
    //////////////////////////////////////
    psi4 = - sol4;
    InsertSlice(psi4, psi, 0   , 0);
    for (int s=1;s<Ls;s++) {
      ExtractSlice(psi4,c,s,0);
       InsertSlice(psi4,psi,s,0);
    }
    
    // Pdag PV^-1 Dm P 
    _Matrix.P(psi,B);
    _Matrix.M(B,A);

    // PauliVillars invert
    _Matrix.SetMass(1.0);    
    _Matrix.Mdag(A,B);    _PauliVillarsSolver(HermOp,B,A);    
    _Matrix.SetMass(m);  

    _Matrix.Pdag(A,B);

    InsertSlice(sol4,B,0,0);

    // Convert from y back to x with P+ circular shift
    _Matrix.P(B,sol5);
    
  }
};

}
}
