#include <Grid.h>
#include <parallelIO/GridNerscIO.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> simd_layout = GridDefaultSimd();
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> latt_size  ({16,16,16,32});
  std::vector<int> clatt_size  ({4,4,4,8});
  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridCartesian     Coarse(clatt_size,simd_layout,mpi_layout);

  LatticeGaugeField Umu(&Fine);

  std::vector<LatticeColourMatrix> U(4,&Fine);
  
  NerscField header;
  
  std::string file("./ckpoint_lat.4000");
  readNerscConfiguration(Umu,header,file);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = peekIndex<3>(Umu,mu);
  }

  // Painful ; fix syntactical niceness
  LatticeComplex LinkTrace(&Fine);
  LinkTrace=zero;
  for(int mu=0;mu<Nd;mu++){
    LinkTrace = LinkTrace + trace(U[mu]);
  }

  // (1+2+3)=6 = N(N-1)/2 terms
  LatticeComplex Plaq(&Fine);
  LatticeComplex cPlaq(&Coarse);
  Plaq = zero;
  for(int mu=1;mu<Nd;mu++){
    for(int nu=0;nu<mu;nu++){
      Plaq = Plaq + trace(U[mu]*Cshift(U[nu],mu,1)*adj(Cshift(U[mu],nu,1))*adj(U[nu]));
    }
  }
  
  double vol = Fine.gSites();
  Complex PlaqScale(1.0/vol/6.0/3.0);
  std::cout <<"PlaqScale" << PlaqScale<<std::endl;

  std::vector<TComplex> Plaq_T(orthosz);
  sliceSum(Plaq,Plaq_T,Nd-1);
  int Nt = Plaq_T.size();

  TComplex Plaq_T_sum; 
  Plaq_T_sum=zero;
  for(int t=0;t<Nt;t++){
    Plaq_T_sum = Plaq_T_sum+Plaq_T[t];
    Complex Pt=TensorRemove(Plaq_T[t]);
    std::cout << "sliced ["<<t<<"]" <<Pt*PlaqScale*Real(Nt)<<std::endl;
  }

  {
    Complex Pt = TensorRemove(Plaq_T_sum);
    std::cout << "total " <<Pt*PlaqScale<<std::endl;
  }  

  TComplex Tp = sum(Plaq);
  Complex p  = TensorRemove(Tp);
  std::cout << "calculated plaquettes " <<p*PlaqScale<<std::endl;


  Complex LinkTraceScale(1.0/vol/4.0/3.0);
  TComplex Tl = sum(LinkTrace);
  Complex l  = TensorRemove(Tl);
  std::cout << "calculated link trace " <<l*LinkTraceScale<<std::endl;

  sumBlocks(cPlaq,Plaq);
  TComplex TcP = sum(cPlaq);
  Complex ll= TensorRemove(TcP);
  std::cout << "coarsened plaquettes sum to " <<ll*PlaqScale<<std::endl;
  
  Grid_finalize();
}
