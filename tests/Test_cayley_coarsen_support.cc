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



  const int Ls=4;

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
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

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
      if (nn!=0) { U[nn]=zero; std::cout << "zeroing gauge field in dir "<<nn<<std::endl; }
      else       { U[nn] = cone;std::cout << "unit gauge field in dir "<<nn<<std::endl; }
    }
    pokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }
#endif  
  RealD mass=0.5;
  RealD M5=1.8;

  DomainWallFermion Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  Gamma5HermitianLinearOperator<DomainWallFermion,LatticeFermion> HermIndefOp(Ddwf);

  HermIndefOp.Op(src,ref);
  HermIndefOp.OpDiag(src,result);
  
  for(int d=0;d<4;d++){
    HermIndefOp.OpDir(src,tmp,d,+1); result=result+tmp; 
    std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
    HermIndefOp.OpDir(src,tmp,d,-1); result=result+tmp;
    std::cout<<"dir "<<d<<" tmp "<<norm2(tmp)<<std::endl;
  }
  err = result-ref;
  std::cout<<"Error "<<norm2(err)<<std::endl;

  const int nbasis = 2;
  std::vector<LatticeFermion> subspace(nbasis,FGrid);
  
  for(int b=0;b<nbasis;b++){
    random(RNG5,subspace[b]);
  }
  std::cout << "Computed randoms"<< std::endl;

  CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOp(*Coarse5d);

  LittleDiracOp.CoarsenOperator(FGrid,HermIndefOp,subspace);

  typedef   Lattice<iVector<vComplex,nbasis > > coarse_vec;
  
  coarse_vec c_src (Coarse5d);  c_src= zero;
  coarse_vec c_res (Coarse5d);
  
  Complex one(1.0);
  c_src = one;  // 1 in every element for vector 1.

  //  TODO
  // -- promote from subspace, check we get the vector we wanted
  // -- apply ldop; check we get the same as inner product of M times big vec
  // -- pick blocks one by one. Evaluate matrix elements.

  std::cout << "Multiplying by LittleDiracOp "<< std::endl;
  LittleDiracOp.M(c_src,c_res);

  std::cout << "Done "<< std::endl;
  Grid_finalize();
}
