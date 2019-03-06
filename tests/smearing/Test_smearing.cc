/*
 *

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/Test_wilson_cg_prec.cc

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


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  LatticeFermion src(&Grid); random(pRNG,src);
  RealD nrm = norm2(src);
  LatticeFermion result(&Grid); result=zero;
  LatticeGaugeField Umu(&Grid);
  //  SU3::HotConfiguration(pRNG,Umu);
  SU3::ColdConfiguration(Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  std::vector<int> site({4,4,0,0});
  src=zero;
  SpinColourVector scv;
  scv=zero;
  scv()(0)(0) = 1.0;
  pokeSite(scv,src,site);

  CovariantSmearing<PeriodicGimplR>::GaussianSmear(U, src, 2.0, 50, Tdir);

  std::cout << src <<std::endl;
  
}

