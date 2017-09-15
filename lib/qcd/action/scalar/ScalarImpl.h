#ifndef SCALAR_IMPL
#define SCALAR_IMPL


namespace Grid {
  //namespace QCD {

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
      GridBase           *grid = out._grid;
      Field              kmu(grid), one(grid);
      const unsigned int nd    = grid->_ndimension;
      std::vector<int>   &l    = grid->_fdimensions;

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
      FFT   fft((GridCartesian *)in._grid);
      Field inFT(in._grid);

      fft.FFT_all_dim(inFT, in, FFT::forward);
      inFT = inFT*momKernel;
      fft.FFT_all_dim(out, inFT, FFT::backward);
    }

    static void FreePropagator(const Field &in, Field &out, RealD m)
    {
      Field momKernel(in._grid);

      MomentumSpacePropagator(momKernel, m);
      FreePropagator(in, out, momKernel);
    }

  };


  #define USE_FFT_ACCELERATION


  template <class S, unsigned int N>
  class ScalarAdjMatrixImplTypes {
  public:
    typedef S Simd;
    typedef QCD::SU<N> Group;

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


    static void MomentaSquare(ComplexField& out){
      GridBase           *grid = out._grid;
      const std::vector<int>   &l    = grid->FullDimensions();
      ComplexField              kmu(grid);
      
      for(int mu = 0; mu < grid->Nd(); mu++)
      {
        Real twoPiL = M_PI*2.0/l[mu];
        LatticeCoordinate(kmu,mu);
        kmu = 2.0*sin(0.5*twoPiL*kmu);
        out += kmu*kmu;
      }
    }

    static void MomentumSpacePropagator(ComplexField &out, RealD m)
    {
      GridBase           *grid = out._grid;
      ComplexField one(grid); one = Complex(1.0,0.0);
      out = m*m;
      MomentaSquare(out);
      out = one/out;
    }


    static inline void generate_momenta(Field& P, GridParallelRNG& pRNG) {
      #ifndef USE_FFT_ACCELERATION
      Group::GaussianFundamentalLieAlgebraMatrix(pRNG, P);
      #else
      
      Field Ptmp(P._grid), Pp(P._grid);
      Group::GaussianFundamentalLieAlgebraMatrix(pRNG, Ptmp);
      // if we change the mass I need a renormalization here
      // transform and multiply by (M*M+p*p)^-1
      GridCartesian *Grid = dynamic_cast<GridCartesian*>(P._grid);
      FFT theFFT(Grid);
      ComplexField p2(Grid);
      RealD M = 1.0;
      p2= zero;

      theFFT.FFT_all_dim(Pp,Ptmp,FFT::forward);
      MomentaSquare(p2);
      p2 += M*M;
      p2 = sqrt(p2);
      Pp *= p2;
      theFFT.FFT_all_dim(P,Pp,FFT::backward);
      
      #endif //USE_FFT_ACCELERATION
    }

    static inline Field projectForce(Field& P) {return P;}

    static inline void update_field(Field& P, Field& U, double ep) {
      #ifndef USE_FFT_ACCELERATION
      U += P*ep;
      #else
      // Here we can eventually add the Fourier acceleration
      // FFT transform P(x) -> P(p)
      // divide by (M^2+p^2)  M external parameter (how to pass?)
      // P'(p) = P(p)/(M^2+p^2)
      // Transform back -> P'(x)
      // U += P'(x)*ep
      
      // the dynamic cast is safe
      GridCartesian *Grid = dynamic_cast<GridCartesian*>(U._grid);
      FFT theFFT(Grid);
      Field Pp(Grid), Pnew(Grid);
      std::vector<int> full_dim = Grid->FullDimensions();

      theFFT.FFT_all_dim(Pp,P,FFT::forward);
      RealD M = 1.0;  
      static bool first_call = true;
      static ComplexField p2(Grid);
      if (first_call){
      MomentumSpacePropagator(p2,M);
      first_call = false;
      }
      Pp *= p2;
      theFFT.FFT_all_dim(Pnew,Pp,FFT::backward);     
      U += Pnew * ep;
      
      #endif //USE_FFT_ACCELERATION
    }

    static inline RealD FieldSquareNorm(Field &U)
    {
      #ifndef USE_FFT_ACCELERATION
      return (TensorRemove(sum(trace(U*U))).real());
      #else
      // In case of Fourier acceleration we have to:
      // compute U(p)*U(p)/(M^2+p^2))   Parseval theorem
      // 1 FFT needed U(x) -> U(p)
      // M to be passed
      
      GridCartesian *Grid = dynamic_cast<GridCartesian *>(U._grid);
      FFT theFFT(Grid);
      Field Up(Grid), Utilde(Grid);
      std::vector<int> full_dim = Grid->FullDimensions();
      
      theFFT.FFT_all_dim(Up, U, FFT::forward);
      RealD M = 1.0;  
      ComplexField p2(Grid);
      MomentumSpacePropagator(p2,M);
      Field Up2 = Up*p2; 
      // from the definition of the DFT we need to divide by the volume
      return (-TensorRemove(sum(trace(adj(Up)*Up2))).real()/U._grid->gSites());
      #endif //USE_FFT_ACCELERATION
    }

    static inline void HotConfiguration(GridParallelRNG &pRNG, Field &U) {
      Group::GaussianFundamentalLieAlgebraMatrix(pRNG, U);
    }

    static inline void TepidConfiguration(GridParallelRNG &pRNG, Field &U) {
      Group::GaussianFundamentalLieAlgebraMatrix(pRNG, U, 0.01);
    }

    static inline void ColdConfiguration(GridParallelRNG &pRNG, Field &U) {
      U = zero;
    }

  };




  typedef ScalarImplTypes<vReal> ScalarImplR;
  typedef ScalarImplTypes<vRealF> ScalarImplF;
  typedef ScalarImplTypes<vRealD> ScalarImplD;
  typedef ScalarImplTypes<vComplex> ScalarImplCR;
  typedef ScalarImplTypes<vComplexF> ScalarImplCF;
  typedef ScalarImplTypes<vComplexD> ScalarImplCD;

  // Hardcoding here the size of the matrices
  typedef ScalarAdjMatrixImplTypes<vComplex,  QCD::Nc> ScalarAdjImplR;
  typedef ScalarAdjMatrixImplTypes<vComplexF, QCD::Nc> ScalarAdjImplF;
  typedef ScalarAdjMatrixImplTypes<vComplexD, QCD::Nc> ScalarAdjImplD;

  template <int Colours > using ScalarNxNAdjImplR = ScalarAdjMatrixImplTypes<vComplex,   Colours >;
  template <int Colours > using ScalarNxNAdjImplF = ScalarAdjMatrixImplTypes<vComplexF,  Colours >;
  template <int Colours > using ScalarNxNAdjImplD = ScalarAdjMatrixImplTypes<vComplexD,  Colours >;

  //}
}

#endif
