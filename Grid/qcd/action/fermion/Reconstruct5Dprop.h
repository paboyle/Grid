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

template<class Field,class PVinverter> class Reconstruct5DfromPhysical {
 private:
  PVinverter & PauliVillarsSolver;
 public:

 /////////////////////////////////////////////////////
 // First cut works, 10 Oct 2018.
 //
 // Must form a plan to get this into production for Zmobius acceleration
 // of the Mobius exact AMA corrections.
 //
 // TODO : understand absence of contact term in eqns in Hantao's thesis
 //        sol4 is contact term subtracted, but thesis & Brower's paper suggests not.
 //
 // Step 1: Localise PV inverse in a routine. [DONE]
 // Step 2: Schur based PV inverse            [DONE]
 // Step 3: Fourier accelerated PV inverse    [DONE]
 //
 /////////////////////////////////////////////////////
 
  Reconstruct5DfromPhysical(PVinverter &_PauliVillarsSolver) 
    : PauliVillarsSolver(_PauliVillarsSolver) 
  { 
  };


   template<class Matrix>
   void PV(Matrix &_Matrix,const Field &src,Field &sol)
   {
     RealD m = _Matrix.Mass();
     _Matrix.SetMass(1.0);
     _Matrix.M(src,sol);
     _Matrix.SetMass(m);
   }
   template<class Matrix>
   void PVdag(Matrix &_Matrix,const Field &src,Field &sol)
   {
     RealD m = _Matrix.Mass();
     _Matrix.SetMass(1.0);
     _Matrix.Mdag(src,sol);
     _Matrix.SetMass(m);
   }
  template<class Matrix>
  void operator() (Matrix & _Matrix,const Field &sol4,const Field &src4, Field &sol5){

    int Ls =  _Matrix.Ls;

    Field psi4(_Matrix.GaugeGrid());
    Field psi(_Matrix.FermionGrid());
    Field A  (_Matrix.FermionGrid());
    Field B  (_Matrix.FermionGrid());
    Field c  (_Matrix.FermionGrid());

    typedef typename Matrix::Coeff_t Coeff_t;

    std::cout << GridLogMessage<< " ************************************************" << std::endl;
    std::cout << GridLogMessage<< " Reconstruct5Dprop: c.f. MADWF algorithm         " << std::endl;
    std::cout << GridLogMessage<< " ************************************************" << std::endl;

    ///////////////////////////////////////
    //Import source, include Dminus factors
    ///////////////////////////////////////
    _Matrix.ImportPhysicalFermionSource(src4,B); 

    ///////////////////////////////////////
    // Set up c from src4
    ///////////////////////////////////////
    PauliVillarsSolver(_Matrix,B,A);
    _Matrix.Pdag(A,c);

    //////////////////////////////////////
    // Build Pdag PV^-1 Dm P [-sol4,c2,c3... cL]
    //////////////////////////////////////
    psi4 = - sol4;
    InsertSlice(psi4, psi, 0   , 0);
    for (int s=1;s<Ls;s++) {
      ExtractSlice(psi4,c,s,0);
       InsertSlice(psi4,psi,s,0);
    }

    /////////////////////////////
    // Pdag PV^-1 Dm P 
    /////////////////////////////
    _Matrix.P(psi,B);
    _Matrix.M(B,A);
    PauliVillarsSolver(_Matrix,A,B);
    _Matrix.Pdag(B,A);

    //////////////////////////////
    // Reinsert surface prop
    //////////////////////////////
    InsertSlice(sol4,A,0,0);

    //////////////////////////////
    // Convert from y back to x 
    //////////////////////////////
    _Matrix.P(A,sol5);
    
  }
};

}
}
