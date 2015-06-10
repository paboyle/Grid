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

double lowpass(double x)
{
  return pow(x*x,-2);
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Chebyshev<LatticeFermion> filter(-120.0, 120.0,256, lowpass);
  ofstream csv(std::string("filter.dat"),std::ios::out|std::ios::trunc);
  filter.csv(csv);
  csv.close();

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

  LatticeFermion    src(FGrid); gaussian(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid); ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid); gaussian(RNG4,Umu);


  //random(RNG4,Umu);
  NerscField header;
  std::string file("./ckpoint_lat.400");
  readNerscConfiguration(Umu,header,file);

#if 0
  LatticeColourMatrix U(UGrid);
  Complex cone(1.0,0.0);
  for(int nn=0;nn<Nd;nn++){
    if (nn==3) { 
      U=zero; std::cout << "zeroing gauge field in dir "<<nn<<std::endl;
    //    else       { U[nn]= cone;std::cout << "unit gauge field in dir "<<nn<<std::endl; }
      pokeIndex<LorentzIndex>(Umu,U,nn);
    }
  }
#endif  
  RealD mass=0.5;
  RealD M5=1.8;

  DomainWallFermion Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  Gamma5R5HermitianLinearOperator<DomainWallFermion,LatticeFermion> HermIndefOp(Ddwf);

  const int nbasis = 8;
  std::vector<LatticeFermion> subspace(nbasis,FGrid);
  LatticeFermion noise(FGrid);
  LatticeFermion ms(FGrid);
  for(int b=0;b<nbasis;b++){
    Gamma g5(Gamma::Gamma5);
    gaussian(RNG5,noise);
    RealD scale = pow(norm2(noise),-0.5); 
    noise=noise*scale;

    HermIndefOp.HermOp(noise,ms); std::cout << "Noise    "<<b<<" Ms "<<norm2(ms)<< " "<< norm2(noise)<<std::endl;

    filter(HermIndefOp,noise,subspace[b]);
    scale = pow(norm2(subspace[b]),-0.5); 
    subspace[b]=subspace[b]*scale;
    HermIndefOp.HermOp(subspace[b],ms); std::cout << "Filtered "<<b<<" Ms "<<norm2(ms)<< " "<<norm2(subspace[b]) <<std::endl;
  }
  std::cout << "Computed randoms"<< std::endl;

  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  LittleDiracOperator LittleDiracOp(*Coarse5d);
  LittleDiracOp.CoarsenOperator(FGrid,HermIndefOp,subspace);
  
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  gaussian(CRNG,c_src);
  c_res=zero;

  std::cout << "Solving MCR on coarse space "<< std::endl;
  HermitianLinearOperator<LittleDiracOperator,CoarseVector> HermIndefLdop(LittleDiracOp);
  ConjugateResidual<CoarseVector> MCR(1.0e-8,10000);
  MCR(HermIndefLdop,c_src,c_res);

  std::cout << "Done "<< std::endl;
  Grid_finalize();
}
