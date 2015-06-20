#include <Grid.h>
#include <qcd/utils/WilsonLoops.h>
#include <qcd/utils/SUn.h>

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
  return pow(x*x+1.0,-2);
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Chebyshev<LatticeFermion> filter(-150.0, 150.0,16, lowpass);
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
  LatticeGaugeField Umu(UGrid); 



  //gaussian(RNG4,Umu);
  //random(RNG4,Umu);

  NerscField header;
  std::string file("./ckpoint_lat.400");
  readNerscConfiguration(Umu,header,file);
  //  SU3::ColdConfiguration(RNG4,Umu);
  //  SU3::TepidConfiguration(RNG4,Umu);
  //  SU3::HotConfiguration(RNG4,Umu);
  //  Umu=zero;

#if 0
  LatticeColourMatrix U(UGrid);
  for(int nn=0;nn<Nd;nn++){
    U=peekIndex<LorentzIndex>(Umu,nn);
    U=U*adj(U)-1.0;
    std::cout<<"SU3 test "<<norm2(U)<<std::endl;
  }
#endif  
  RealD mass=0.1;
  RealD M5=1.5;

  DomainWallFermion Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  Gamma5R5HermitianLinearOperator<DomainWallFermion,LatticeFermion> HermIndefOp(Ddwf);

  const int nbasis = 8;

#if 0
  std::vector<LatticeFermion> subspace(nbasis,FGrid);
  LatticeFermion noise(FGrid);
  LatticeFermion ms(FGrid);
  for(int b=0;b<nbasis;b++){

    gaussian(RNG5,noise);
    RealD scale = pow(norm2(noise),-0.5); 
    noise=noise*scale;

    HermIndefOp.Op(noise,ms); std::cout << "Noise    "<<b<<" Ms "<<norm2(ms)<< " "<< norm2(noise)<<std::endl;

    //    filter(HermIndefOp,noise,subspace[b]);
    // inverse iteration
    MdagMLinearOperator<DomainWallFermion,LatticeFermion> HermDefOp(Ddwf);
    ConjugateGradient<LatticeFermion> CG(1.0e-4,10000);

    for(int i=0;i<1;i++){

      CG(HermDefOp,noise,subspace[b]);
      noise = subspace[b];

      scale = pow(norm2(noise),-0.5); 
      noise=noise*scale;
      HermDefOp.Op(noise,ms); std::cout << "filt    "<<b<<" <u|H|u> "<<norm2(ms)<< " "<< norm2(noise)<<std::endl;
    }

    subspace[b] = noise;
    HermIndefOp.Op(subspace[b],ms); std::cout << "Filtered "<<b<<" Ms "<<norm2(ms)<< " "<<norm2(subspace[b]) <<std::endl;

  }
  std::cout << "Computed randoms"<< std::endl;
#else
  std::cout<<"Calling Aggregation class" <<std::endl;
  MdagMLinearOperator<DomainWallFermion,LatticeFermion> HermDefOp(Ddwf);
  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FGrid);
  Aggregates.CreateSubspace(RNG5,HermDefOp);
  std::cout << "Called aggregation class"<< std::endl;
#endif


  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  LittleDiracOperator LittleDiracOp(*Coarse5d);
  LittleDiracOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);
  
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  gaussian(CRNG,c_src);
  c_res=zero;

  std::cout << "Solving CG on coarse space "<< std::endl;
  MdagMLinearOperator<LittleDiracOperator,CoarseVector> PosdefLdop(LittleDiracOp);
  ConjugateGradient<CoarseVector> CG(1.0e-6,10000);
  CG(PosdefLdop,c_src,c_res);

  std::cout << "Done "<< std::endl;
  Grid_finalize();
}
