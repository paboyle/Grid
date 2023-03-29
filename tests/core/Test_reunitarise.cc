    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_quenched_update.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt({8,8,8,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(latt, 
							GridDefaultSimd(Nd,vComplexD::Nsimd()),
							GridDefaultMpi());

  GridCartesian * gridF = SpaceTimeGrid::makeFourDimGrid(latt, 
							GridDefaultSimd(Nd,vComplexF::Nsimd()),
							GridDefaultMpi());
  

  ///////////////////////////////
  // Configuration of known size
  ///////////////////////////////
  LatticeColourMatrixD ident(grid);
  LatticeColourMatrixD U(grid);
  LatticeColourMatrixD UU(grid);
  LatticeColourMatrixD tmp(grid);
  LatticeColourMatrixD org(grid);
  LatticeColourMatrixF UF(gridF);

  LatticeGaugeField Umu(grid);

  ident =1.0;

  // RNG set up for test
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  SU<Nc>::HotConfiguration(pRNG,Umu);

  U = PeekIndex<LorentzIndex>(Umu,0);
  org=U;


  tmp=  U*adj(U) - ident ;
  RealD Def1 = norm2( tmp );
  std::cout << " Defect1 "<<Def1<<std::endl;

  tmp = U - org;
  std::cout << "Diff1 "<<norm2(tmp)<<std::endl;
  precisionChange(UF,U);
  precisionChange(U,UF);

  tmp=  U*adj(U) - ident ;
  RealD Def2 = norm2(  tmp );
  std::cout << " Defect2 "<<Def2<<std::endl;

  tmp = U - org;
  std::cout << "Diff2 "<<norm2(tmp)<<std::endl;

  U = ProjectOnGroup(U);

  tmp=  U*adj(U) - ident ;
  RealD Def3 = norm2(  tmp);
  std::cout << " Defect3 "<<Def3<<std::endl;


  tmp = U - org;
  std::cout << "Diff3 "<<norm2(tmp)<<std::endl;

  LatticeComplexD detU(grid);
  LatticeComplexD detUU(grid);

  detU= Determinant(U) ;
  detU=detU-1.0;
  std::cout << "Determinant defect before screw up " << norm2(detU)<<std::endl;

  std::cout << " Screwing up determinant " << std::endl;

  RealD theta = 0.2;
  ComplexD phase(cos(theta),sin(theta));
  for(int i=0;i<Nc;i++){
    auto element = PeekIndex<ColourIndex>(U,Nc-1,i);
    element = element * phase;
    PokeIndex<ColourIndex>(U,element,Nc-1,i);
  }
  U=U*0.1;
  UU=U;

  detU= Determinant(U) ;
  detU=detU-1.0;
  std::cout << "Determinant defect before projection " <<norm2(detU)<<std::endl;
  tmp = U*adj(U) - ident;
  std::cout << "Unitarity check before projection    " << norm2(tmp)<<std::endl; 
#if (Nc == 3)
  ProjectSU3(U);
  detU= Determinant(U) ;
  detU= detU -1.0;
  std::cout << "Determinant ProjectSU3 defect " <<norm2(detU)<<std::endl;
  tmp = U*adj(U) - ident;
  std::cout << "Unitarity check after projection    " << norm2(tmp)<<std::endl; 
#endif
  
  ProjectSUn(UU);
  detUU= Determinant(UU);
  detUU= detUU -1.0;
  std::cout << "Determinant ProjectSUn defect " <<norm2(detUU)<<std::endl;
  tmp = UU*adj(UU) - ident;
  std::cout << "Unitarity check after projection    " << norm2(tmp)<<std::endl; 
  
  Grid_finalize();
}




