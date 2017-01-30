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
using namespace Grid::QCD;

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

int toint(const char* str){
  std::stringstream os; os << str;
  int out; os >> out;
  return out;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  assert(argc >= 5);
  
  std::vector<int> latt(4,0);
  latt[0] = toint(argv[1]);
  latt[1] = toint(argv[2]);
  latt[2] = toint(argv[3]);
  latt[3] = toint(argv[4]);
  
  const int Ls= toint(argv[5]);
  
  std::cout << "Lattice size (" << latt[0] << "," << latt[1] << "," << latt[2] << "," << latt[3] << ") Ls=" << Ls << std::endl;
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
  std::cout << "SIMD layout (" << simd_layout[0] << "," << simd_layout[1] << "," << simd_layout[2] << "," << simd_layout[3] << ")" << std::endl;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt, simd_layout,GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  //typedef Lattice<iGparitySpinColourVector<vComplexD> > LatticeType;
  typedef LatticeFermionD LatticeType;
  
  LatticeType    src(FGrid); random(RNG5,src);

  LatticeType src_o(FrbGrid);
  pickCheckerboard(Odd,src_o,src);

  std::vector<int> site(5);
  std::vector<int> cbsite(5);
  typedef typename GridTypeMapper<LatticeType::vector_object>::scalar_object sobj;

  // std::cout << "sizeof(vobj) " << sizeof(LatticeType::vector_object) << std::endl;
  // std::cout << "sizeof(sobj) " << sizeof(sobj) << std::endl;
  std::cout << "v1 from uncheckerboarded field, v2 from odd-parity red-black field\n";

  for(site[0]=0;site[0]<Ls;site[0]++){
    for(site[4]=0;site[4]<latt[3];site[4]++){
      for(site[3]=0;site[3]<latt[2];site[3]++){
	for(site[2]=0;site[2]<latt[1];site[2]++){
	  for(site[1]=0;site[1]<latt[0];site[1]++){
	    if(src_o._grid->CheckerBoard(site) != src_o.checkerboard)
	      continue;

	    std::cout << "Site (" << site[0] << "," << site[1] << "," << site[2] << "," << site[3] << "," << site[4] << ")" << std::endl;
	    sobj v1, v2;
	    peekLocalSite(v1,src,site);
	    peekLocalSite(v2,src_o,site);
	    
	    RealD v1_norm = norm2(v1);
	    RealD v2_norm = norm2(v2);
	    RealD diff = v2_norm - v1_norm;

	    std::cout << v1_norm << " " << v2_norm << " " << diff << '\n';
	    if(fabs(diff) > 1e-12){
	      std::cout << "ERROR!\n";
	      exit(-1);
	    }
	  }
	}
      }
    }
  }

  





  
  // LatticeFermion result(FGrid); result=zero;
  // LatticeGaugeField Umu(UGrid); 

  // SU3::HotConfiguration(RNG4,Umu);

  // std::vector<LatticeColourMatrix> U(4,UGrid);
  // for(int mu=0;mu<Nd;mu++){
  //   U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  // }
  
  // RealD mass=0.1;
  // RealD M5=1.8;
  // DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  // LatticeFermion    src_o(FrbGrid);
  // LatticeFermion result_o(FrbGrid);
  // pickCheckerboard(Odd,src_o,src);
  // result_o=zero;

  // SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOpEO(Ddwf);
  // ConjugateGradient<LatticeFermion> CG(1.0e-8,10000);
  // CG(HermOpEO,src_o,result_o);

  Grid_finalize();
}
