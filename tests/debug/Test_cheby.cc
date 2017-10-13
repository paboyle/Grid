    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cheby.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

RealD InverseApproximation(RealD x){
  return 1.0/x;
}
RealD SqrtApproximation(RealD x){
  return std::sqrt(x);
}
RealD Approximation32(RealD x){
  return std::pow(x,-1.0/32.0);
}
RealD Approximation2(RealD x){
  return std::pow(x,-1.0/2.0);
}

RealD StepFunction(RealD x){
  if ( x<10.0 )  return 1.0;
  else return 0.0;
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
						       GridDefaultSimd(Nd,vComplex::Nsimd()),
						       GridDefaultMpi());

  double     lo=0.1;
  double     hi=64.0;

  Chebyshev<LatticeFermion> ChebyInv(lo,hi,2000,InverseApproximation);

  {
    std::ofstream of("chebyinv");
    ChebyInv.csv(of);
  }

  ChebyInv.JacksonSmooth();

  {
    std::ofstream of("chebyinvjack");
    ChebyInv.csv(of);
  }


  Chebyshev<LatticeFermion> ChebyStep(lo,hi,200,StepFunction);
  {
    std::ofstream of("chebystep");
    ChebyStep.csv(of);
  }


  ChebyStep.JacksonSmooth();
  {
    std::ofstream of("chebystepjack");
    ChebyStep.csv(of);
  }

  lo=-8;
  hi=8;
  Chebyshev<LatticeFermion> ChebyIndefInv(lo,hi,40,InverseApproximation);
  {
    std::ofstream of("chebyindefinv");
    ChebyIndefInv.csv(of);
  }

  lo=0;
  hi=64;
  Chebyshev<LatticeFermion> ChebyNE(lo,hi,40,InverseApproximation);
  {
    std::ofstream of("chebyNE");
    ChebyNE.csv(of);
  }

  lo=0.0;
  hi=4.0;
  Chebyshev<LatticeFermion> Cheby32(lo,hi,2000,Approximation32);
  {
    std::ofstream of("cheby32");
    Cheby32.csv(of);
  }
  Cheby32.JacksonSmooth();
  {
    std::ofstream of("cheby32jack");
    Cheby32.csv(of);
  }

  Chebyshev<LatticeFermion> ChebySqrt(lo,hi,2000,Approximation2);
  {
    std::ofstream of("chebysqrt");
    ChebySqrt.csv(of);
  }
  ChebySqrt.JacksonSmooth();
  {
    std::ofstream of("chebysqrtjack");
    ChebySqrt.csv(of);
  }


  Grid_finalize();
}
