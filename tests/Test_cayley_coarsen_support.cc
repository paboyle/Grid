    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_cayley_coarsen_support.cc

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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


template<class d>
struct scal {
  d internal;
};

  Gamma::GammaMatrix Gmu [] = {
    Gamma::GammaX,
    Gamma::GammaY,
    Gamma::GammaZ,
    Gamma::GammaT
  };

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid
  std::vector<int> clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/2;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  LatticeFermion    src(FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid); ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);

#if 0
  std::vector<LatticeColourMatrix> U(4,UGrid);
  Umu=zero;
  Complex cone(1.0,0.0);
  for(int nn=0;nn<Nd;nn++){
    if(1) {
      if (nn>2) { U[nn]=zero; std::cout<<GridLogMessage << "zeroing gauge field in dir "<<nn<<std::endl; }
      else      { U[nn]=cone; std::cout<<GridLogMessage << "unit gauge field in dir "<<nn<<std::endl; }
    }
    pokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
#endif  

  RealD mass=0.5;
  RealD M5=1.8;

  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> HermIndefOp(Ddwf);

  HermIndefOp.Op(src,ref);
  HermIndefOp.OpDiag(src,result);
  
  for(int d=0;d<4;d++){
    HermIndefOp.OpDir(src,tmp,d+1,+1); result=result+tmp; 
    std::cout<<GridLogMessage<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    HermIndefOp.OpDir(src,tmp,d+1,-1); result=result+tmp;
    std::cout<<GridLogMessage<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
  }
  err = result-ref;
  std::cout<<GridLogMessage<<"Error "<<norm2(err)<<std::endl;

  const int nbasis = 2;
  LatticeFermion prom(FGrid);

  std::vector<LatticeFermion> subspace(nbasis,FGrid);

  std::cout<<GridLogMessage<<"Calling Aggregation class" <<std::endl;

  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermDefOp(Ddwf);
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FGrid);
  Aggregates.CreateSubspaceRandom(RNG5);

  subspace=Aggregates.subspace;

  std::cout<<GridLogMessage << "Called aggregation class"<< std::endl;

  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  LittleDiracOperator LittleDiracOp(*Coarse5d);

  LittleDiracOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);
  
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  CoarseVector c_proj(Coarse5d);
  
  Complex one(1.0);
  c_src = one;  // 1 in every element for vector 1.
  
  blockPromote(c_src,err,subspace);

  prom=zero;
  for(int b=0;b<nbasis;b++){
    prom=prom+subspace[b];
  }
  err=err-prom; 
  std::cout<<GridLogMessage<<"Promoted back from subspace err "<<norm2(err)<<std::endl;

  HermIndefOp.HermOp(prom,tmp);
  blockProject(c_proj,tmp,subspace);

  LittleDiracOp.M(c_src,c_res);

  c_proj = c_proj - c_res;
  std::cout<<GridLogMessage<<"Representation of ldop within subspace "<<norm2(c_proj)<<std::endl;

  std::cout<<GridLogMessage << "Multiplying by LittleDiracOp "<< std::endl;
  LittleDiracOp.M(c_src,c_res);

  std::cout<<GridLogMessage<<"Testing hermiticity explicitly by inspecting matrix elements"<<std::endl;
  LittleDiracOp.AssertHermitian();

  std::cout<<GridLogMessage << "Testing Hermiticity stochastically "<< std::endl;
  CoarseVector phi(Coarse5d);
  CoarseVector chi(Coarse5d);
  CoarseVector Aphi(Coarse5d);
  CoarseVector Achi(Coarse5d);

  random(CRNG,phi);
  random(CRNG,chi);


  std::cout<<GridLogMessage<<"Made randoms"<<std::endl;

  LittleDiracOp.M(phi,Aphi);
  LittleDiracOp.Mdag(chi,Achi);

  ComplexD pAc = innerProduct(chi,Aphi);
  ComplexD cAp = innerProduct(phi,Achi);
  ComplexD cAc = innerProduct(chi,Achi);
  ComplexD pAp = innerProduct(phi,Aphi);

  std::cout<<GridLogMessage<< "pAc "<<pAc<<" cAp "<< cAp<< " diff "<<pAc-adj(cAp)<<std::endl;

  std::cout<<GridLogMessage<< "pAp "<<pAp<<" cAc "<< cAc<<"Should be real"<< std::endl;

  std::cout<<GridLogMessage<<"Testing linearity"<<std::endl;
  CoarseVector PhiPlusChi(Coarse5d);
  CoarseVector APhiPlusChi(Coarse5d);
  CoarseVector linerr(Coarse5d);
  PhiPlusChi = phi+chi;
  LittleDiracOp.M(PhiPlusChi,APhiPlusChi);

  linerr= APhiPlusChi-Aphi;
  linerr= linerr-Achi;
  std::cout<<GridLogMessage<<"**Diff "<<norm2(linerr)<<std::endl;


  std::cout<<GridLogMessage << "Done "<< std::endl;
  Grid_finalize();
}
