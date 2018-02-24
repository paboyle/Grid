#ifndef SCALAR_IMPL
#define SCALAR_IMPL

NAMESPACE_BEGIN(Grid);

template <class S>
class ScalarImplTypes {
public:
  typedef S Simd;

  template <typename vtype>
  using iImplField = iScalar<iScalar<iScalar<vtype> > >;

  typedef iImplField<Simd> SiteField;
  typedef SiteField        SitePropagator;
  typedef SiteField        SiteComplex;
    
  typedef Lattice<SiteField> Field;
  typedef Field              ComplexField;
  typedef Field              FermionField;
  typedef Field              PropagatorField;
    
  static inline void generate_momenta(Field& P, GridParallelRNG& pRNG){
    gaussian(pRNG, P);
  }

  static inline Field projectForce(Field& P){return P;}

  static inline void update_field(Field& P, Field& U, double ep) {
    U += P*ep;
  }

  static inline RealD FieldSquareNorm(Field& U) {
    return (- sum(trace(U*U))/2.0);
  }

  static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
    gaussian(pRNG, U);
  }

  static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
    gaussian(pRNG, U);
  }

  static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
    U = 1.0;
  }
    
  static void MomentumSpacePropagator(Field &out, RealD m)
  {
    GridBase           *grid = out.Grid();
    Field              kmu(grid), one(grid);
    const unsigned int nd    = grid->_ndimension;
    Coordinate         &l    = grid->_fdimensions;
      
    one = Complex(1.0,0.0);
    out = m*m;
    for(int mu = 0; mu < nd; mu++)
      {
        Real twoPiL = M_PI*2./l[mu];
        
        LatticeCoordinate(kmu,mu);
        kmu = 2.*sin(.5*twoPiL*kmu);
        out = out + kmu*kmu;
      }
    out = one/out;
  }
    
  static void FreePropagator(const Field &in, Field &out,
			     const Field &momKernel)
  {
    FFT   fft((GridCartesian *)in.Grid());
    Field inFT(in.Grid());
      
    fft.FFT_all_dim(inFT, in, FFT::forward);
    inFT = inFT*momKernel;
    fft.FFT_all_dim(out, inFT, FFT::backward);
  }
    
  static void FreePropagator(const Field &in, Field &out, RealD m)
  {
    Field momKernel(in.Grid());
      
    MomentumSpacePropagator(momKernel, m);
    FreePropagator(in, out, momKernel);
  }
    
};

template <class S, unsigned int N>
class ScalarAdjMatrixImplTypes {
public:
  typedef S Simd;
  typedef SU<N> Group;
    
  template <typename vtype>
  using iImplField   = iScalar<iScalar<iMatrix<vtype, N>>>;
  template <typename vtype>
  using iImplComplex = iScalar<iScalar<iScalar<vtype>>>;

  typedef iImplField<Simd>   SiteField;
  typedef SiteField          SitePropagator;
  typedef iImplComplex<Simd> SiteComplex;
    
  typedef Lattice<SiteField>   Field;
  typedef Lattice<SiteComplex> ComplexField;
  typedef Field                FermionField;
  typedef Field                PropagatorField;

  static inline void generate_momenta(Field& P, GridParallelRNG& pRNG) {
    Group::GaussianFundamentalLieAlgebraMatrix(pRNG, P);
  }

  static inline Field projectForce(Field& P) {return P;}

  static inline void update_field(Field& P, Field& U, double ep) {
    U += P*ep;
  }

  static inline RealD FieldSquareNorm(Field& U) {
    return (TensorRemove(sum(trace(U*U))).real());
  }

  static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
    Group::GaussianFundamentalLieAlgebraMatrix(pRNG, U);
  }

  static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
    Group::GaussianFundamentalLieAlgebraMatrix(pRNG, U, 0.01);
  }

  static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
    U = Zero();
  }

};




typedef ScalarImplTypes<vReal> ScalarImplR;
typedef ScalarImplTypes<vRealF> ScalarImplF;
typedef ScalarImplTypes<vRealD> ScalarImplD;
typedef ScalarImplTypes<vComplex> ScalarImplCR;
typedef ScalarImplTypes<vComplexF> ScalarImplCF;
typedef ScalarImplTypes<vComplexD> ScalarImplCD;
    
// Hardcoding here the size of the matrices
typedef ScalarAdjMatrixImplTypes<vComplex,  Nc> ScalarAdjImplR;
typedef ScalarAdjMatrixImplTypes<vComplexF, Nc> ScalarAdjImplF;
typedef ScalarAdjMatrixImplTypes<vComplexD, Nc> ScalarAdjImplD;

template <int Colours > using ScalarNxNAdjImplR = ScalarAdjMatrixImplTypes<vComplex,   Colours >;
template <int Colours > using ScalarNxNAdjImplF = ScalarAdjMatrixImplTypes<vComplexF,  Colours >;
template <int Colours > using ScalarNxNAdjImplD = ScalarAdjMatrixImplTypes<vComplexD,  Colours >;

NAMESPACE_END(Grid);

#endif
