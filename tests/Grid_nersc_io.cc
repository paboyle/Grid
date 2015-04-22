#include <Grid.h>
#include <parallelIO/GridNerscIO.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> simd_layout({1,1,2,2});
  std::vector<int> mpi_layout ({2,1,1,2});
  std::vector<int> latt_size  ({16,16,16,32});
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridRNG           FineRNG(&Fine);
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
  Plaq = zero;
  for(int mu=1;mu<Nd;mu++){
    for(int nu=0;nu<mu;nu++){
      Plaq = Plaq + trace(U[mu]*Cshift(U[nu],mu,1)*adj(Cshift(U[mu],nu,1))*adj(U[nu]));
    }
  }


  double vol = Fine.gSites();

  Complex PlaqScale(1.0/vol/6.0/3.0);
  TComplex Tp = sum(Plaq);
  Complex p  = TensorRemove(Tp);
  std::cout << "calculated plaquettes " <<p*PlaqScale<<std::endl;

  Complex LinkTraceScale(1.0/vol/4.0/3.0);
  TComplex Tl = sum(LinkTrace);
  Complex l  = TensorRemove(Tl);
  std::cout << "calculated link trace " <<l*LinkTraceScale<<std::endl;

  Grid_finalize();
}
