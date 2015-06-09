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

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid
  std::vector<int> clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/2;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());;
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
      if (nn>2) { U[nn]=zero; std::cout << "zeroing gauge field in dir "<<nn<<std::endl; }
      else      { U[nn]=cone; std::cout << "unit gauge field in dir "<<nn<<std::endl; }
    }
    pokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
#endif  

  RealD mass=0.5;
  RealD M5=1.8;

  DomainWallFermion Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  Gamma5R5HermitianLinearOperator<DomainWallFermion,LatticeFermion> HermIndefOp(Ddwf);

  HermIndefOp.Op(src,ref);
  HermIndefOp.OpDiag(src,result);
  
  for(int d=0;d<4;d++){
    HermIndefOp.OpDir(src,tmp,d+1,+1); result=result+tmp; 
    std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    HermIndefOp.OpDir(src,tmp,d+1,-1); result=result+tmp;
    std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
  }
  err = result-ref;
  std::cout<<"Error "<<norm2(err)<<std::endl;

  const int nbasis = 2;
  std::vector<LatticeFermion> subspace(nbasis,FGrid);
  LatticeFermion prom(FGrid);

  for(int b=0;b<nbasis;b++){
    random(RNG5,subspace[b]);
  }
  std::cout << "Computed randoms"<< std::endl;

  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  LittleDiracOperator LittleDiracOp(*Coarse5d);

  LittleDiracOp.CoarsenOperator(FGrid,HermIndefOp,subspace);
  
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  CoarseVector c_proj(Coarse5d);
  
  //  TODO
  // -- promote from subspace, check we get the vector we wanted
  // -- apply ldop; check we get the same as inner product of M times big vec
  // -- pick blocks one by one. Evaluate matrix elements.
  Complex one(1.0);
  c_src = one;  // 1 in every element for vector 1.
  
  blockPromote(c_src,err,subspace);


  prom=zero;
  for(int b=0;b<nbasis;b++){
    prom=prom+subspace[b];
  }
  err=err-prom; 
  std::cout<<"Promoted back from subspace err "<<norm2(err)<<std::endl;

  HermIndefOp.HermOp(prom,tmp);
  blockProject(c_proj,tmp,subspace);

  LittleDiracOp.M(c_src,c_res);

  c_proj = c_proj - c_res;
  std::cout<<"Representation of ldop within subspace "<<norm2(c_proj)<<std::endl;

  std::cout << "Multiplying by LittleDiracOp "<< std::endl;
  LittleDiracOp.M(c_src,c_res);

  std::cout<<"Testing hermiticity explicitly by inspecting matrix elements"<<std::endl;
  LittleDiracOp.AssertHermitian();

  std::cout << "Testing Hermiticity stochastically "<< std::endl;
  CoarseVector phi(Coarse5d);
  CoarseVector chi(Coarse5d);
  CoarseVector Aphi(Coarse5d);
  CoarseVector Achi(Coarse5d);

  random(CRNG,phi);
  random(CRNG,chi);


  std::cout<<"Made randoms"<<std::endl;

  LittleDiracOp.M(phi,Aphi);
  LittleDiracOp.Mdag(chi,Achi);

  ComplexD pAc = innerProduct(chi,Aphi);
  ComplexD cAp = innerProduct(phi,Achi);
  ComplexD cAc = innerProduct(chi,Achi);
  ComplexD pAp = innerProduct(phi,Aphi);

  std::cout<< "pAc "<<pAc<<" cAp "<< cAp<< " diff "<<pAc-adj(cAp)<<std::endl;

  std::cout<< "pAp "<<pAp<<" cAc "<< cAc<<"Should be real"<< std::endl;

  std::cout<<"Testing linearity"<<std::endl;
  CoarseVector PhiPlusChi(Coarse5d);
  CoarseVector APhiPlusChi(Coarse5d);
  CoarseVector linerr(Coarse5d);
  PhiPlusChi = phi+chi;
  LittleDiracOp.M(PhiPlusChi,APhiPlusChi);

  linerr= APhiPlusChi-Aphi;
  linerr= linerr-Achi;
  std::cout<<"**Diff "<<norm2(linerr)<<std::endl;


  std::cout << "Done "<< std::endl;
  Grid_finalize();
}
