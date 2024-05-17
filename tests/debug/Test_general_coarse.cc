    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_padded_cell.cc

    Copyright (C) 2023

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
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

using namespace std;
using namespace Grid;

gridblasHandle_t GridBLAS::gridblasHandle;
int            GridBLAS::gridblasInit;

///////////////////////
// Tells little dirac op to use MdagM as the .Op()
///////////////////////
template<class Field>
class HermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> & wrapped;
public:
  HermOpAdaptor(LinearOperatorBase<Field> &wrapme) : wrapped(wrapme)  {};
  void OpDiag (const Field &in, Field &out) {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){    assert(0);  };
  void Op     (const Field &in, Field &out){
    wrapped.HermOp(in,out);
  }
  void AdjOp     (const Field &in, Field &out){
    wrapped.HermOp(in,out);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
  void HermOp(const Field &in, Field &out){
    wrapped.HermOp(in,out);
  }
};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=4;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid
  Coordinate clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/4;
  }

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt,
							    GridDefaultSimd(Nd,vComplex::Nsimd()),
							    GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  LatticeFermion    src(FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid); ref=Zero();
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  //  Umu=Zero();
  
  RealD mass=0.1;
  RealD M5=1.8;

  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  const int nbasis = 62;
  const int cb = 0 ;
  LatticeFermion prom(FGrid);

  std::vector<LatticeFermion> subspace(nbasis,FGrid);

  std::cout<<GridLogMessage<<"Calling Aggregation class" <<std::endl;

  ///////////////////////////////////////////////////////////
  // Squared operator is in HermOp
  ///////////////////////////////////////////////////////////
  MdagMLinearOperator<DomainWallFermionD,LatticeFermion> HermDefOp(Ddwf);

  ///////////////////////////////////////////////////
  // Random aggregation space
  ///////////////////////////////////////////////////
  std::cout<<GridLogMessage << "Building random aggregation class"<< std::endl;
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FGrid,cb);
  Aggregates.CreateSubspaceRandom(RNG5);

  ///////////////////////////////////////////////////
  // Build little dirac op
  ///////////////////////////////////////////////////
  std::cout<<GridLogMessage << "Building little Dirac operator"<< std::endl;

  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);
  LittleDiracOperator LittleDiracOp(geom,FGrid,Coarse5d);
  LittleDiracOperator LittleDiracOpCol(geom,FGrid,Coarse5d);

  HermOpAdaptor<LatticeFermionD> HOA(HermDefOp);

  LittleDiracOp.CoarsenOperator(HOA,Aggregates);
  
  ///////////////////////////////////////////////////
  // Test the operator
  ///////////////////////////////////////////////////
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  CoarseVector c_res_dag(Coarse5d);
  CoarseVector c_proj(Coarse5d);

  subspace=Aggregates.subspace;

  //  random(CRNG,c_src);
  c_src = 1.0;

  blockPromote(c_src,err,subspace);

  prom=Zero();
  for(int b=0;b<nbasis;b++){
    prom=prom+subspace[b];
  }
  err=err-prom; 
  std::cout<<GridLogMessage<<"Promoted back from subspace: err "<<norm2(err)<<std::endl;
  std::cout<<GridLogMessage<<"c_src "<<norm2(c_src)<<std::endl;
  std::cout<<GridLogMessage<<"prom  "<<norm2(prom)<<std::endl;

  HermDefOp.HermOp(prom,tmp);

  blockProject(c_proj,tmp,subspace);
  std::cout<<GridLogMessage<<" Called Big Dirac Op "<<norm2(tmp)<<std::endl;

  std::cout<<GridLogMessage<<" Calling little Dirac Op "<<std::endl;
  LittleDiracOp.M(c_src,c_res);
  LittleDiracOp.Mdag(c_src,c_res_dag);

  std::cout<<GridLogMessage<<"Little dop : "<<norm2(c_res)<<std::endl;
  std::cout<<GridLogMessage<<"Little dop dag : "<<norm2(c_res_dag)<<std::endl;
  std::cout<<GridLogMessage<<"Big dop in subspace : "<<norm2(c_proj)<<std::endl;

  c_proj = c_proj - c_res;
  std::cout<<GridLogMessage<<" ldop error: "<<norm2(c_proj)<<std::endl;

  c_res_dag = c_res_dag - c_res;
  std::cout<<GridLogMessage<<"Little dopDag - dop: "<<norm2(c_res_dag)<<std::endl;

  std::cout<<GridLogMessage << "Testing Hermiticity stochastically "<< std::endl;
  CoarseVector phi(Coarse5d);
  CoarseVector chi(Coarse5d);
  CoarseVector Aphi(Coarse5d);
  CoarseVector Achi(Coarse5d);

  random(CRNG,phi);
  random(CRNG,chi);

  std::cout<<GridLogMessage<<"Made randoms "<<norm2(phi)<<" " << norm2(chi)<<std::endl;

  LittleDiracOp.M(phi,Aphi);

  LittleDiracOp.Mdag(chi,Achi);

  std::cout<<GridLogMessage<<"Aphi "<<norm2(Aphi)<<" A chi" << norm2(Achi)<<std::endl;

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

  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////////
  //  Create a higher dim coarse grid
  //////////////////////////////////////////////////////////////////////////////////////

  const int nrhs=vComplex::Nsimd()*3;

  Coordinate mpi=GridDefaultMpi();
  Coordinate rhMpi ({1,1,mpi[0],mpi[1],mpi[2],mpi[3]});
  Coordinate rhLatt({nrhs,1,clatt[0],clatt[1],clatt[2],clatt[3]});
  Coordinate rhSimd({vComplex::Nsimd(),1, 1,1,1,1});

  GridCartesian *CoarseMrhs = new GridCartesian(rhLatt,rhSimd,rhMpi); 

  
  MultiGeneralCoarsenedMatrix mrhs(LittleDiracOp,CoarseMrhs);
  typedef decltype(mrhs) MultiGeneralCoarsenedMatrix_t;
  
  //////////////////////////////////////////
  // Test against single RHS
  //////////////////////////////////////////
  {
    GridParallelRNG          rh_CRNG(CoarseMrhs);rh_CRNG.SeedFixedIntegers(cseeds);
    CoarseVector rh_phi(CoarseMrhs);
    CoarseVector rh_res(CoarseMrhs);
    random(rh_CRNG,rh_phi);

    std::cout << "Warmup"<<std::endl;
    mrhs.M(rh_phi,rh_res);
    const int ncall=5;
    RealD t0=-usecond();
    for(int i=0;i<ncall;i++){
      std::cout << "Call "<<i<<"/"<<ncall<<std::endl;
      mrhs.M(rh_phi,rh_res);
    }
    t0+=usecond();
    RealD t1=0;
    for(int r=0;r<nrhs;r++){
      std::cout << " compare to single RHS "<<r<<"/"<<nrhs<<std::endl;
      ExtractSlice(phi,rh_phi,r,0);
      ExtractSlice(chi,rh_res,r,0);
      LittleDiracOp.M(phi,Aphi);
      t1-=usecond();
      for(int i=0;i<ncall;i++){
	std::cout << "Call "<<i<<"/"<<ncall<<std::endl;
	LittleDiracOp.M(phi,Aphi);
      }
      t1+=usecond();
      Coordinate site({0,0,0,0,0});
      auto  bad = peekSite(chi,site);
      auto good = peekSite(Aphi,site);
      std::cout << " mrhs [" <<r <<"] "<< norm2(chi)<<std::endl;
      std::cout << " srhs [" <<r <<"] "<< norm2(Aphi)<<std::endl;
      chi=chi-Aphi;
      RealD diff =norm2(chi);
      std::cout << r << " diff " << diff<<std::endl;
      assert(diff < 1.0e-10);
    }
    std::cout << nrhs<< " mrhs " << t0/ncall/nrhs <<" us"<<std::endl;
    std::cout << nrhs<< " srhs " << t1/ncall/nrhs <<" us"<<std::endl;
  }

  //////////////////////////////////////////
  // Test against single RHS
  //////////////////////////////////////////
  {
    typedef HermitianLinearOperator<MultiGeneralCoarsenedMatrix_t,CoarseVector> HermMatrix;
    HermMatrix MrhsCoarseOp     (mrhs);

    GridParallelRNG          rh_CRNG(CoarseMrhs);rh_CRNG.SeedFixedIntegers(cseeds);
    ConjugateGradient<CoarseVector>  mrhsCG(1.0e-8,2000,true);
    CoarseVector rh_res(CoarseMrhs);
    CoarseVector rh_src(CoarseMrhs);
    random(rh_CRNG,rh_src);
    rh_res= Zero();
    mrhsCG(MrhsCoarseOp,rh_src,rh_res);
  }
  
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;
  std::cout<<GridLogMessage<<"*******************************************"<<std::endl;

  Grid_finalize();
  return 0;
}
