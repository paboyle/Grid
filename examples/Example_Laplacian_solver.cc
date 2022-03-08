#include <Grid/Grid.h>
using namespace Grid;

template<class Field>
void SimpleConjugateGradient(LinearOperatorBase<Field> &HPDop,const Field &b, Field &x)
{
    RealD cp, c, alpha, d, beta, ssq;
    RealD Tolerance=1.0e-10;
    int MaxIterations=10000;
    
    Field p(b), mmp(b), r(b);

    HPDop.HermOpAndNorm(x, mmp, d, beta);
    
    r = b - mmp;
    p = r;

    cp = alpha = norm2(p);
    ssq = norm2(b);

    RealD rsq = Tolerance * Tolerance * ssq;

    for (int k = 1; k <= MaxIterations; k++) {
      c = cp;

      HPDop.HermOp(p, mmp);

      d = real(innerProduct(p,mmp));

      alpha = c / d;

      r = r - alpha *mmp;
      cp = norm2(r);
      beta = cp / c;

      x   = x   + alpha* p ;
      p   = r   + beta* p ;

      std::cout << "iteration "<<k<<" cp " <<std::sqrt(cp/ssq) << std::endl;
      if (cp <= rsq) {
        return;
      }
    }
    assert(0);
}



// Flip sign to make prop to p^2, not -p^2 relative to last example
template<class Gimpl,class Field> class CovariantLaplacianCshift : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  GridBase *grid;
  GaugeField U;
  RealD m2=1.0e-2;
  CovariantLaplacianCshift(GaugeField &_U)    :
    grid(_U.Grid()),
    U(_U) {  };

  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &in, Field &out)
  {
    out=Zero();
    for(int mu=0;mu<Nd-1;mu++) {
      GaugeLinkField Umu = PeekIndex<LorentzIndex>(U, mu); // NB: Inefficent
      out = out - Gimpl::CovShiftForward(Umu,mu,in);    
      out = out - Gimpl::CovShiftBackward(Umu,mu,in);    
      out = out + 2.0*in + m2*in;
    }
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};



int main(int argc, char ** argv)
{
  Grid_init(&argc, &argv);

  typedef LatticeColourVector Field;

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();

  GridCartesian    Grid(latt_size,simd_layout,mpi_layout);
  GridParallelRNG  RNG(&Grid);  RNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));


  LatticeGaugeField U(&Grid);

  SU<Nc>::ColdConfiguration(RNG,U);

  typedef CovariantLaplacianCshift <PeriodicGimplR,Field> Laplacian_t;
  Laplacian_t Laplacian(U);

  
  ColourVector ColourKronecker;
  ColourKronecker = Zero();
  ColourKronecker()()(0) = 1.0;

  Coordinate site({0,0,0,0}); // Point source at origin

  Field kronecker(&Grid);
  kronecker = Zero();
  pokeSite(ColourKronecker,kronecker,site);

  Field psi(&Grid); psi=Zero();

  HermitianLinearOperator<Laplacian_t,Field> HermOp(Laplacian);
  SimpleConjugateGradient(HermOp, kronecker,psi);

  Field r(&Grid);
  Laplacian.M(psi,r);
  r=kronecker-r;
  
  std::cout << "True residual "<< norm2(r) <<std::endl;

  // Optionally print the result vector
  // std::cout << psi<<std::endl;

  Grid_finalize();
}
