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


template<class What> 
void  TestWhat(What & Ddwf,
	       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
	       GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
	       RealD mass, RealD M5,
	       GridParallelRNG *RNG4,   GridParallelRNG *RNG5);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  const int Ls=9;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);


  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  RealD mass=0.1;
  RealD M5  =1.8;
  std::cout <<"ContinuedFractionFermion test"<<std::endl;
  ContinuedFractionFermion5D Dcf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  TestWhat<ContinuedFractionFermion5D>(Dcf,FGrid,FrbGrid,UGrid,UrbGrid,mass,M5,&RNG4,&RNG5);

  Grid_finalize();
}

template<class What> 
void  TestWhat(What & Ddwf, 
	       GridCartesian         * FGrid,	       GridRedBlackCartesian * FrbGrid,
	       GridCartesian         * UGrid,	       GridRedBlackCartesian * UrbGrid,
	       RealD mass, RealD M5,
	       GridParallelRNG *RNG4,
	       GridParallelRNG *RNG5)
{

  LatticeFermion src   (FGrid); random(*RNG5,src);
  LatticeFermion phi   (FGrid); random(*RNG5,phi);
  LatticeFermion chi   (FGrid); random(*RNG5,chi);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid);    ref=zero;
  LatticeFermion    tmp(FGrid);    tmp=zero;
  LatticeFermion    err(FGrid);    tmp=zero;

  LatticeFermion src_e (FrbGrid);
  LatticeFermion src_o (FrbGrid);
  LatticeFermion r_e   (FrbGrid);
  LatticeFermion r_o   (FrbGrid);
  LatticeFermion r_eo  (FGrid);
  LatticeFermion r_eeoo(FGrid);

  std::cout<<"=========================================================="<<std::endl;
  std::cout<<"= Testing that Meo + Moe + Moo + Mee = Munprec "<<std::endl;
  std::cout<<"=========================================================="<<std::endl;

  pickCheckerboard(Even,src_e,src);
  pickCheckerboard( Odd,src_o,src);

  Ddwf.Meooe(src_e,r_o);  std::cout<<"Applied Meo "<<norm2(r_o)<<std::endl;
  Ddwf.Meooe(src_o,r_e);  std::cout<<"Applied Moe "<<norm2(r_e)<<std::endl;
  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  Ddwf.Mooee(src_e,r_e);  std::cout<<"Applied Mee"<<norm2(r_e)<<std::endl;
  Ddwf.Mooee(src_o,r_o);  std::cout<<"Applied Moo"<<norm2(r_o)<<std::endl;
  setCheckerboard(r_eeoo,r_e);
  setCheckerboard(r_eeoo,r_o);

  r_eo=r_eo+r_eeoo;
  Ddwf.M(src,ref);  

  //  std::cout << r_eo<<std::endl;
  //  std::cout << ref <<std::endl;

  err= ref - r_eo;
  std::cout << "EO norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(r_eo) <<std::endl;
    
  LatticeComplex cerr(FGrid);
  cerr = localInnerProduct(err,err);
  //  std::cout << cerr<<std::endl;

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test Ddagger is the dagger of D by requiring                "<<std::endl;
  std::cout<<"=  < phi | Deo | chi > * = < chi | Deo^dag| phi>  "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;
  
  LatticeFermion chi_e   (FrbGrid);
  LatticeFermion chi_o   (FrbGrid);

  LatticeFermion dchi_e  (FrbGrid);
  LatticeFermion dchi_o  (FrbGrid);

  LatticeFermion phi_e   (FrbGrid);
  LatticeFermion phi_o   (FrbGrid);

  LatticeFermion dphi_e  (FrbGrid);
  LatticeFermion dphi_o  (FrbGrid);


  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);

  Ddwf.Meooe(chi_e,dchi_o);
  Ddwf.Meooe(chi_o,dchi_e);
  Ddwf.MeooeDag(phi_e,dphi_o);
  Ddwf.MeooeDag(phi_o,dphi_e);

  ComplexD pDce = innerProduct(phi_e,dchi_e);
  ComplexD pDco = innerProduct(phi_o,dchi_o);
  ComplexD cDpe = innerProduct(chi_e,dphi_e);
  ComplexD cDpo = innerProduct(chi_o,dphi_o);

  std::cout <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout <<"pDce - conj(cDpo) "<< pDce-conj(cDpo) <<std::endl;
  std::cout <<"pDco - conj(cDpe) "<< pDco-conj(cDpe) <<std::endl;

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test MeeInv Mee = 1                                         "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ddwf.Mooee(chi_e,src_e);
  Ddwf.MooeeInv(src_e,phi_e);

  Ddwf.Mooee(chi_o,src_o);
  Ddwf.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Ddwf.MooeeDag(chi_e,src_e);
  Ddwf.MooeeInvDag(src_e,phi_e);

  Ddwf.MooeeDag(chi_o,src_o);
  Ddwf.MooeeInvDag(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test MpcDagMpc is Hermitian              "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;
  
  random(*RNG5,phi);
  random(*RNG5,chi);
  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);
  RealD t1,t2;

  Ddwf.MpcDagMpc(chi_e,dchi_e,t1,t2);
  Ddwf.MpcDagMpc(chi_o,dchi_o,t1,t2);

  Ddwf.MpcDagMpc(phi_e,dphi_e,t1,t2);
  Ddwf.MpcDagMpc(phi_o,dphi_o,t1,t2);

  pDce = innerProduct(phi_e,dchi_e);
  pDco = innerProduct(phi_o,dchi_o);
  cDpe = innerProduct(chi_e,dphi_e);
  cDpo = innerProduct(chi_o,dphi_o);

  std::cout <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout <<"pDce - conj(cDpo) "<< pDco-conj(cDpo) <<std::endl;
  std::cout <<"pDco - conj(cDpe) "<< pDce-conj(cDpe) <<std::endl;
  
}
