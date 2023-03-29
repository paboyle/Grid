/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/core/Test_prec_change.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Ls = 12;
  Coordinate latt4 = GridDefaultLatt();

  GridCartesian         * UGridD   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  GridCartesian         * FGridD   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridD);
  GridRedBlackCartesian * FrbGridD = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridD);

  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);

  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  
  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGridD);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG5F(FGridF);  RNG5F.SeedFixedIntegers(seeds5);
  std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

  LatticeFermionD field_d(FGridD), tmp_d(FGridD);
  random(RNG5,field_d);
  RealD norm2_field_d = norm2(field_d);
  
  LatticeFermionD2 field_d2(FGridF), tmp_d2(FGridF);
  random(RNG5F,field_d2);
  RealD norm2_field_d2 = norm2(field_d2);
  
  LatticeFermionF field_f(FGridF);
  
  //Test original implementation
  {
    std::cout << GridLogMessage << "Testing original implementation" << std::endl;
    field_f = Zero();
    precisionChangeOrig(field_f,field_d);
    RealD Ndiff = (norm2_field_d - norm2(field_f))/norm2_field_d;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of single and double prec fields differs by " << Ndiff << std::endl;
    tmp_d = Zero();
    precisionChangeOrig(tmp_d, field_f);
    Ndiff = norm2( LatticeFermionD(tmp_d-field_d) ) / norm2_field_d;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of back-converted and original double prec fields differs by " << Ndiff << std::endl;
  }
  //Test new implementation with pregenerated workspace
  {
    std::cout << GridLogMessage << "Testing new implementation with pregenerated workspace" << std::endl;
    precisionChangeWorkspace wk_sp_to_dp(field_d.Grid(),field_f.Grid());
    precisionChangeWorkspace wk_dp_to_sp(field_f.Grid(),field_d.Grid());
    
    field_f = Zero();
    precisionChange(field_f,field_d,wk_dp_to_sp);
    RealD Ndiff = (norm2_field_d - norm2(field_f))/norm2_field_d;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of single and double prec fields differs by " << Ndiff << std::endl;
    tmp_d = Zero();
    precisionChange(tmp_d, field_f,wk_sp_to_dp);
    Ndiff = norm2( LatticeFermionD(tmp_d-field_d) ) / norm2_field_d;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of back-converted and original double prec fields differs by " << Ndiff << std::endl;
  }
  //Test new implementation without pregenerated workspace
  {
    std::cout << GridLogMessage << "Testing new implementation without pregenerated workspace" << std::endl;
    field_f = Zero();
    precisionChange(field_f,field_d);
    RealD Ndiff = (norm2_field_d - norm2(field_f))/norm2_field_d;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of single and double prec fields differs by " << Ndiff << std::endl;
    tmp_d = Zero();
    precisionChange(tmp_d, field_f);
    Ndiff = norm2( LatticeFermionD(tmp_d-field_d) ) / norm2_field_d;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of back-converted and original double prec fields differs by " << Ndiff << std::endl;
  } 
  //Test fast implementation
  {
    std::cout << GridLogMessage << "Testing fast (double2) implementation" << std::endl;
    field_f = Zero();
    precisionChangeFast(field_f,field_d2);
    RealD Ndiff = (norm2_field_d2 - norm2(field_f))/norm2_field_d2;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of single and double prec fields differs by " << Ndiff << std::endl;
    tmp_d2 = Zero();
    precisionChangeFast(tmp_d2, field_f);
    Ndiff = norm2( LatticeFermionD2(tmp_d2-field_d2) ) / norm2_field_d2;
    std::cout << GridLogMessage << (fabs(Ndiff) > 1e-05 ? "!!FAIL" : "Pass") << ": relative norm2 of back-converted and original double prec fields differs by " << Ndiff << std::endl;
  }
  std::cout << "Done" << std::endl;
  
  Grid_finalize();
}
