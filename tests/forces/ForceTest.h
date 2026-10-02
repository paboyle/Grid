/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/forces/ForceTest.h

    Copyright (C) 2022

Author: Peter Boyle <pboyle@bnl.gov>

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

/////////////////////////////////////////////////////////////////////////////////////////////
// Finite-difference check of an action's force: refresh, S1 = S(U), step the links by eps
// along a random momentum P, take the force at the midpoint, step again, S2 = S(U'').
// The force predicts S2 - S1 to O(eps^3).
//
// Runs on a copy of U, so the caller's field is unchanged. Returns S2 - S1 - dSpred for
// the caller to assert on; halving eps should reduce it by about 8.
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Gimpl>
RealD ForceTest(Action<LatticeGaugeField> &action,LatticeGaugeField & Uin,MomentumFilterBase<LatticeGaugeField> &Filter,RealD eps=0.005)
{
  GridBase *UGrid = Uin.Grid();

  std::vector<int> seeds({1,2,3,5});
  GridSerialRNG            sRNG;         sRNG.SeedFixedIntegers(seeds);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds);

  LatticeColourMatrix Pmu(UGrid);
  LatticeGaugeField P(UGrid);
  LatticeGaugeField UdSdU(UGrid);
  LatticeGaugeField U(UGrid);
  U = Uin;

  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
  std::cout << GridLogMessage << " Force test for "<<action.action_name()<<" eps "<<eps<<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Refresh "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  Gimpl::generate_momenta(P,sRNG,RNG4);
  Filter.applyFilter(P);

  action.refresh(U,sRNG,RNG4);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  RealD S1 = action.S(U);

  Gimpl::update_field(P,U,eps);

  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Derivative "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  action.deriv(U,UdSdU);
  UdSdU = Ta(UdSdU);
  Filter.applyFilter(UdSdU);

  DumpSliceNorm("Force",UdSdU,Nd-1);

  Gimpl::update_field(P,U,eps);
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout << GridLogMessage << " Action "<<action.action_name()<<std::endl;
  std::cout << GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  RealD S2 = action.S(U);

  // Use the derivative
  LatticeComplex dS(UGrid); dS = Zero();
  for(int mu=0;mu<Nd;mu++){
    auto UdSdUmu = PeekIndex<LorentzIndex>(UdSdU,mu);
    Pmu= PeekIndex<LorentzIndex>(P,mu);
    dS = dS - trace(Pmu*UdSdUmu)*eps*2.0*2.0;
  }
  ComplexD dSpred    = sum(dS);
  RealD diff =  S2-S1-dSpred.real();

  std::cout<< GridLogMessage << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<< GridLogMessage << "S1 : "<< S1    <<std::endl;
  std::cout<< GridLogMessage << "S2 : "<< S2    <<std::endl;
  std::cout<< GridLogMessage << "dS : "<< S2-S1 <<std::endl;
  std::cout<< GridLogMessage << "dSpred : "<< dSpred.real() <<std::endl;
  std::cout<< GridLogMessage << "diff : "<< diff<<std::endl;
  std::cout<< GridLogMessage << "*********************************************************"<<std::endl;
  std::cout<< GridLogMessage << "Done" <<std::endl;
  std::cout << GridLogMessage << "*********************************************************"<<std::endl;
  return diff;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// ForceTest at eps0, eps0/2, ... (neps values). A correct force leaves an O(eps^3)
// discrepancy, so successive ratios approach 8 as eps shrinks; an error in the force adds
// an O(eps) part and drives them towards 2. A strongly curved action needs a smaller eps0
// before the ratios settle. Returns the ratio for the smallest pair.
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Gimpl>
RealD ForceTestScaling(Action<LatticeGaugeField> &action,LatticeGaugeField &U,MomentumFilterBase<LatticeGaugeField> &Filter,const std::string &name,
                       RealD eps0=0.01,int neps=4)
{
  std::vector<RealD> eps(neps);
  std::vector<RealD> d(neps);
  for(int i=0;i<neps;i++){
    eps[i] = eps0/(1<<i);
    d[i]   = ForceTest<Gimpl>(action,U,Filter,eps[i]);
  }
  for(int i=0;i<neps;i++){
    std::cout << GridLogMessage << name << ": eps " << eps[i] << " discrepancy " << d[i];
    if ( i > 0 ) {
      std::cout << "  ratio to previous " << d[i-1]/d[i];
    }
    std::cout << std::endl;
  }
  RealD scale = d[neps-2]/d[neps-1];
  std::cout << GridLogMessage << name << ": smallest-eps ratio " << scale << " (expect ~8)" << std::endl;
  return scale;
}
