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

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexF::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(latt_size,simd_layout,mpi_layout);

  int threads = GridThread::GetThreads();
  std::cout << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  //  std::vector<int> seeds({1,2,3,4});
  //  pRNG.SeedFixedIntegers(seeds);
  pRNG.SeedRandomDevice();

  LatticeFermion src   (&Grid); random(pRNG,src);
  LatticeFermion phi   (&Grid); random(pRNG,phi);
  LatticeFermion chi   (&Grid); random(pRNG,chi);
  LatticeFermion result(&Grid); result=zero;
  LatticeFermion    ref(&Grid);    ref=zero;
  LatticeFermion    tmp(&Grid);    tmp=zero;
  LatticeFermion    err(&Grid);    tmp=zero;
  LatticeGaugeField Umu(&Grid); random(pRNG,Umu);
  std::vector<LatticeColourMatrix> U(4,&Grid);

  double volume=1;
  for(int mu=0;mu<Nd;mu++){
    volume=volume*latt_size[mu];
  }  

  // Only one non-zero (y)
  Umu=zero;
  for(int nn=0;nn<Nd;nn++){
    random(pRNG,U[nn]);
    pokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }

  RealD mass=0.1;

  WilsonFermion Dw(Umu,Grid,RBGrid,mass);

  LatticeFermion src_e   (&RBGrid);
  LatticeFermion src_o   (&RBGrid);
  LatticeFermion r_e   (&RBGrid);
  LatticeFermion r_o   (&RBGrid);
  LatticeFermion r_eo  (&Grid);

  std::cout<<"=========================================================="<<std::endl;
  std::cout<<"= Testing that Deo + Doe = Dunprec "<<std::endl;
  std::cout<<"=========================================================="<<std::endl;

  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  Dw.Meooe(src_e,r_o);  std::cout<<"Applied Meo"<<std::endl;
  Dw.Meooe(src_o,r_e);  std::cout<<"Applied Moe"<<std::endl;
  Dw.Dhop (src,ref,DaggerNo);

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err= ref - r_eo;
  std::cout << "EO norm diff   "<< norm2(err)<< " "<<norm2(ref)<< " " << norm2(r_eo) <<std::endl;

  LatticeComplex cerr(&Grid);
  cerr = localInnerProduct(err,err);

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test Ddagger is the dagger of D by requiring                "<<std::endl;
  std::cout<<"=  < phi | Deo | chi > * = < chi | Deo^dag| phi>  "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;
  
  LatticeFermion chi_e   (&RBGrid);
  LatticeFermion chi_o   (&RBGrid);

  LatticeFermion dchi_e  (&RBGrid);
  LatticeFermion dchi_o  (&RBGrid);

  LatticeFermion phi_e   (&RBGrid);
  LatticeFermion phi_o   (&RBGrid);

  LatticeFermion dphi_e  (&RBGrid);
  LatticeFermion dphi_o  (&RBGrid);


  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);

  Dw.Meooe(chi_e,dchi_o);
  Dw.Meooe(chi_o,dchi_e);
  Dw.MeooeDag(phi_e,dphi_o);
  Dw.MeooeDag(phi_o,dphi_e);

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

  Dw.Mooee(chi_e,src_e);
  Dw.MooeeInv(src_e,phi_e);

  Dw.Mooee(chi_o,src_o);
  Dw.MooeeInv(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test MeeInvDag MeeDag = 1                                   "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;

  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);

  Dw.MooeeDag(chi_e,src_e);
  Dw.MooeeInvDag(src_e,phi_e);

  Dw.MooeeDag(chi_o,src_o);
  Dw.MooeeInvDag(src_o,phi_o);
  
  setCheckerboard(phi,phi_e);
  setCheckerboard(phi,phi_o);

  err = phi-chi;
  std::cout << "norm diff   "<< norm2(err)<< std::endl;

  std::cout<<"=============================================================="<<std::endl;
  std::cout<<"= Test MpcDagMpc is Hermitian              "<<std::endl;
  std::cout<<"=============================================================="<<std::endl;
  
  random(pRNG,phi);
  random(pRNG,chi);
  pickCheckerboard(Even,chi_e,chi);
  pickCheckerboard(Odd ,chi_o,chi);
  pickCheckerboard(Even,phi_e,phi);
  pickCheckerboard(Odd ,phi_o,phi);
  RealD t1,t2;

  Dw.MpcDagMpc(chi_e,dchi_e,t1,t2);
  Dw.MpcDagMpc(chi_o,dchi_o,t1,t2);

  Dw.MpcDagMpc(phi_e,dphi_e,t1,t2);
  Dw.MpcDagMpc(phi_o,dphi_o,t1,t2);

  pDce = innerProduct(phi_e,dchi_e);
  pDco = innerProduct(phi_o,dchi_o);
  cDpe = innerProduct(chi_e,dphi_e);
  cDpo = innerProduct(chi_o,dphi_o);

  std::cout <<"e "<<pDce<<" "<<cDpe <<std::endl;
  std::cout <<"o "<<pDco<<" "<<cDpo <<std::endl;

  std::cout <<"pDce - conj(cDpo) "<< pDco-conj(cDpo) <<std::endl;
  std::cout <<"pDco - conj(cDpe) "<< pDce-conj(cDpe) <<std::endl;
  
  Grid_finalize();
}
