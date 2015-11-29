#include <fenv.h>
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

static int
FEenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
    old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}


template<class Field> class DumbOperator  : public LinearOperatorBase<Field> {
public:
  LatticeComplex scale;

  DumbOperator(GridBase *grid)    : scale(grid)
  {
    GridParallelRNG  pRNG(grid);  
    std::vector<int> seeds({5,6,7,8});
    pRNG.SeedFixedIntegers(seeds);

    random(pRNG,scale);

    scale = exp(-real(scale)*3.0);
    std::cout << " True matrix \n"<< scale <<std::endl;
  }

  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {};
  void OpDir  (const Field &in, Field &out,int dir,int disp){};

  void Op     (const Field &in, Field &out){
    out = scale * in;
  }
  void AdjOp  (const Field &in, Field &out){
    out = scale * in;
  }
  void HermOp(const Field &in, Field &out){
    double n1, n2;
    HermOpAndNorm(in,out,n1,n2);
  }
  void HermOpAndNorm(const Field &in, Field &out,double &n1,double &n2){
    ComplexD dot;

    out = scale * in;

    dot= innerProduct(in,out);
    n1=real(dot);

    dot = innerProduct(out,out);
    n2=real(dot);
  }
};


int main (int argc, char ** argv)
{

  FEenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 

  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
						       GridDefaultSimd(Nd,vComplex::Nsimd()),
						       GridDefaultMpi());

  GridParallelRNG  RNG(grid);  
  std::vector<int> seeds({1,2,3,4});
  RNG.SeedFixedIntegers(seeds);


  RealD alpha = 1.0;
  RealD beta  = 0.03;
  RealD mu    = 0.0;
  int order = 11;
  ChebyshevLanczos<LatticeComplex> Cheby(alpha,beta,mu,order);
  std::ofstream file("cheby.dat");
  Cheby.csv(file);

  HermOpOperatorFunction<LatticeComplex> X;
  DumbOperator<LatticeComplex> HermOp(grid);

  const int Nk = 40;
  const int Nm = 80;
  const int Nit= 10000;

  int Nconv;
  RealD eresid = 1.0e-8;

  ImplicitlyRestartedLanczos<LatticeComplex> IRL(HermOp,X,Nk,Nm,eresid,Nit);

  ImplicitlyRestartedLanczos<LatticeComplex> ChebyIRL(HermOp,Cheby,Nk,Nm,eresid,Nit);

  LatticeComplex src(grid); gaussian(RNG,src);
  {
    std::vector<RealD>          eval(Nm);
    std::vector<LatticeComplex> evec(Nm,grid);
    IRL.calc(eval,evec,src, Nconv);
  }
  
  {
    //    std::vector<RealD>          eval(Nm);
    //    std::vector<LatticeComplex> evec(Nm,grid);
    //    ChebyIRL.calc(eval,evec,src, Nconv);
  }

  Grid_finalize();
}
