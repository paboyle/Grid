    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/core/Test_precision_change.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>

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

using namespace Grid;


int main (int argc, char ** argv){
  Grid_init(&argc, &argv);
  int Ls = 16;
  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << " and Ls=" << Ls << std::endl;
  GridCartesian* UGrid_d = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridCartesian* FGrid_d = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid_d);
  GridRedBlackCartesian* FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_d);

  GridCartesian* UGrid_f = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridCartesian* FGrid_f = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid_f);
  GridRedBlackCartesian* FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_f);


  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid_d);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid_d);
  RNG4.SeedFixedIntegers(seeds4);

  //Gauge fields
  LatticeGaugeFieldD Umu_d(UGrid_d);
  LatticeGaugeFieldF Umu_f(UGrid_f);
  LatticeGaugeFieldD Umu_d_r(UGrid_d);
  LatticeGaugeFieldD Utmp_d(UGrid_d);

  for(int i=0;i<5;i++){
    random(RNG4, Umu_d);

    precisionChange(Umu_f, Umu_d);
    std::cout << GridLogMessage << "Norm of double-prec and single-prec gauge fields (should be ~equal): " << norm2(Umu_d) << " " << norm2(Umu_f) << std::endl;
    precisionChange(Umu_d_r, Umu_f);
    RealD normdiff = axpy_norm(Utmp_d, -1.0, Umu_d_r, Umu_d);
    std::cout << GridLogMessage << "Norm of difference of back-converted double-prec gauge fields (should be ~0) = " << normdiff << std::endl;
  }

  //Fermion fields
  LatticeFermionD psi_d(FGrid_d);
  LatticeFermionF psi_f(FGrid_f);
  LatticeFermionD psi_d_r(FGrid_d);
  LatticeFermionD psi_tmp_d(FGrid_d);

  for(int i=0;i<5;i++){
    random(RNG5, psi_d);

    precisionChange(psi_f, psi_d);
    std::cout << GridLogMessage << "Norm of double-prec and single-prec fermion fields (should be ~equal): " << norm2(psi_d) << " " << norm2(psi_f) << std::endl;
    precisionChange(psi_d_r, psi_f);
    RealD normdiff = axpy_norm(psi_tmp_d, -1.0, psi_d_r, psi_d);
    std::cout << GridLogMessage << "Norm of difference of back-converted double-prec fermion fields (should be ~0)= " << normdiff << std::endl;
  }

  //Checkerboarded fermion fields
  LatticeFermionD psi_cb_d(FrbGrid_d);
  LatticeFermionF psi_cb_f(FrbGrid_f);
  LatticeFermionD psi_cb_d_r(FrbGrid_d);
  LatticeFermionD psi_cb_tmp_d(FrbGrid_d);

  for(int i=0;i<5;i++){
    random(RNG5, psi_d);
    pickCheckerboard(Odd, psi_cb_d, psi_d);
     
    precisionChange(psi_cb_f, psi_cb_d);
    std::cout << GridLogMessage << "Norm of odd-cb double-prec and single-prec fermion fields (should be ~equal): " << norm2(psi_cb_d) << " " << norm2(psi_cb_f) << std::endl;
    precisionChange(psi_cb_d_r, psi_cb_f);
    RealD normdiff = axpy_norm(psi_cb_tmp_d, -1.0, psi_cb_d_r, psi_cb_d);
    std::cout << GridLogMessage << "Norm of difference of back-converted odd-cb double-prec fermion fields (should be ~0)= " << normdiff << std::endl;


    pickCheckerboard(Even, psi_cb_d, psi_d);
     
    precisionChange(psi_cb_f, psi_cb_d);
    std::cout << GridLogMessage << "Norm of even-cb double-prec and single-prec fermion fields (should be ~equal): " << norm2(psi_cb_d) << " " << norm2(psi_cb_f) << std::endl;
    precisionChange(psi_cb_d_r, psi_cb_f);
    normdiff = axpy_norm(psi_cb_tmp_d, -1.0, psi_cb_d_r, psi_cb_d);
    std::cout << GridLogMessage << "Norm of difference of back-converted even-cb double-prec fermion fields (should be ~0)= " << normdiff << std::endl;
  }



  Grid_finalize();
}
