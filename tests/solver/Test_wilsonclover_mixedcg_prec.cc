    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

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
 ;

template<class d>
struct scal {
  d internal;
};

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout << GridLogMessage << "::::: NB: to enable a quick bit reproducibility check use the --checksums flag. " << std::endl;

  GridCartesian         *FGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridCartesian         *FGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid_d);
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid_f);

  std::vector<int> fSeeds({1, 2, 3, 4});
  GridParallelRNG  fPRNG(FGrid_d);
  fPRNG.SeedFixedIntegers(fSeeds);

  // clang-format off
  LatticeFermionD    src(FGrid_d);    gaussian(fPRNG, src);
  LatticeFermionD    result(FGrid_d); result = Zero();
  LatticeGaugeFieldD Umu_d(FGrid_d);  SU3::HotConfiguration(fPRNG, Umu_d);
  LatticeGaugeFieldF Umu_f(FGrid_f);  precisionChange(Umu_f, Umu_d);
  // clang-format on
  
  RealD mass = -0.1;
  RealD csw_r = 1.0;
  RealD csw_t = 1.0;
  WilsonCloverFermionD Dw_d(Umu_d, *FGrid_d, *FrbGrid_d, mass, csw_r, csw_t);
  WilsonCloverFermionF Dw_f(Umu_f, *FGrid_f, *FrbGrid_f, mass, csw_r, csw_t);

  LatticeFermionD      src_o(FrbGrid_d);
  LatticeFermionD   result_o(FrbGrid_d);
  LatticeFermionD result_o_2(FrbGrid_d);
  pickCheckerboard(Odd, src_o, src);
  result_o.Checkerboard() = Odd;
  result_o = Zero();
  result_o_2.Checkerboard() = Odd;
  result_o_2 = Zero();

  SchurDiagMooeeOperator<WilsonCloverFermionD, LatticeFermionD> HermOpEO_d(Dw_d);
  SchurDiagMooeeOperator<WilsonCloverFermionF, LatticeFermionF> HermOpEO_f(Dw_f);

  std::cout << GridLogMessage << "::::::::::::: Starting mixed CG" << std::endl;
  MixedPrecisionConjugateGradient<LatticeFermionD, LatticeFermionF> mCG(1.0e-8, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO_d);
  mCG(src_o, result_o);

  std::cout << GridLogMessage << "::::::::::::: Starting regular CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8, 10000);
  CG(HermOpEO_d, src_o, result_o_2);

  LatticeFermionD diff_o(FrbGrid_d);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << GridLogMessage << "::::::::::::: Diff between mixed and regular CG: " << diff << std::endl;

  #ifdef HAVE_LIME
  if( GridCmdOptionExists(argv,argv+argc,"--checksums") )
  {
    std::string file1("./Propagator1");
    emptyUserRecord record;
    uint32_t nersc_csum;
    uint32_t scidac_csuma;
    uint32_t scidac_csumb;
    typedef SpinColourVectorD   FermionD;
    typedef vSpinColourVectorD vFermionD;

    BinarySimpleMunger<FermionD,FermionD> munge;
    std::string format = getFormatString<vFermionD>();
    
    BinaryIO::writeLatticeObject<vFermionD,FermionD>(result_o,file1,munge, 0, format,
                 nersc_csum,scidac_csuma,scidac_csumb);

    std::cout << GridLogMessage << " Mixed checksums "<<std::hex << scidac_csuma << " "<<scidac_csumb<<std::endl;

    BinaryIO::writeLatticeObject<vFermionD,FermionD>(result_o_2,file1,munge, 0, format,
                 nersc_csum,scidac_csuma,scidac_csumb);

    std::cout << GridLogMessage << " CG checksums "<<std::hex << scidac_csuma << " "<<scidac_csumb<<std::endl;
  }
  #endif
  
  Grid_finalize();
}
