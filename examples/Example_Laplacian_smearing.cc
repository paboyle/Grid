#include <Grid/Grid.h>
using namespace Grid;

// Function used for Chebyshev smearing
// 
Real MomentumSmearing(Real p2)
{
  return (1 - 4.0*p2) * exp(-p2/4);
}
Real DistillationSmearing(Real p2)
{
  if ( p2 > 0.5 ) return 0.0;
  else return 1.0;
}

// Flip sign to make prop to p^2, not -p^2 relative to last example
template<class Gimpl,class Field> class CovariantLaplacianCshift : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  GridBase *grid;
  GaugeField U;
  
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
      out = out + 2.0*in;
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

  Coordinate site({latt_size[0]/2,
                   latt_size[1]/2,
                   latt_size[2]/2,
		   0});

  Field kronecker(&Grid);
  kronecker = Zero();
  pokeSite(ColourKronecker,kronecker,site);


  Field psi(&Grid), chi(&Grid);

  //////////////////////////////////////
  // Classic Wuppertal smearing
  //////////////////////////////////////

  Integer Iterations = 80;
  Real width = 2.0;
  Real coeff = (width*width) / Real(4*Iterations);

  chi=kronecker;
  //  chi = (1-p^2/2N)^N kronecker
  for(int n = 0; n < Iterations; ++n) {
    Laplacian.M(chi,psi);
    chi = chi - coeff*psi;
  }

  std::cout << " Wuppertal smeared operator is chi = \n" << chi <<std::endl;

  /////////////////////////////////////
  // Chebyshev smearing
  /////////////////////////////////////
  RealD lo = 0.0;
  RealD hi = 12.0; // Analytic free field bound
  HermitianLinearOperator<Laplacian_t,Field> HermOp(Laplacian);

  std::cout << " Checking spectral range of our POSITIVE definite operator \n";
  PowerMethod<Field> PM;
  PM(HermOp,kronecker);
  
  //  Chebyshev<Field> ChebySmear(lo,hi,20,DistillationSmearing);
  Chebyshev<Field> ChebySmear(lo,hi,20,MomentumSmearing);
  {
    std::ofstream of("chebysmear");
    ChebySmear.csv(of);
  }

  ChebySmear(HermOp,kronecker,chi);
  
  std::cout << " Chebyshev smeared operator is chi = \n" << chi <<std::endl;

  Grid_finalize();
}
