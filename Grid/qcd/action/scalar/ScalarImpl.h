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

  #ifdef  USE_FFT_ACCELERATION
  #ifndef FFT_MASS
  #error  "USE_FFT_ACCELERATION is defined but not FFT_MASS"
  #endif
  #endif
  
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

    static void MomentaSquare(ComplexField &out)
    {
      GridBase *grid = out._grid;
      const std::vector<int> &l = grid->FullDimensions();
      ComplexField kmu(grid);

      for (int mu = 0; mu < grid->Nd(); mu++)
      {
        Real twoPiL = M_PI * 2.0 / l[mu];
        LatticeCoordinate(kmu, mu);
        kmu = 2.0 * sin(0.5 * twoPiL * kmu);
        out += kmu * kmu;
      }
    }

    static void MomentumSpacePropagator(ComplexField &out, RealD m)
    {
      GridBase *grid = out._grid;
      ComplexField one(grid);
      one = Complex(1.0, 0.0);
      out = m * m;
      MomentaSquare(out);
      out = one / out;
    }

    static inline void generate_momenta(Field &P, GridParallelRNG &pRNG)
    {
#ifndef USE_FFT_ACCELERATION
      Group::GaussianFundamentalLieAlgebraMatrix(pRNG, P);
#else

      Field Pgaussian(P._grid), Pp(P._grid);
      ComplexField p2(P._grid); p2 = zero;
      RealD M = FFT_MASS;
      
      Group::GaussianFundamentalLieAlgebraMatrix(pRNG, Pgaussian);

      FFT theFFT((GridCartesian*)P._grid);
      theFFT.FFT_all_dim(Pp, Pgaussian, FFT::forward);
      MomentaSquare(p2);
      p2 += M * M;
      p2 = sqrt(p2);
      Pp *= p2;
      theFFT.FFT_all_dim(P, Pp, FFT::backward);

#endif //USE_FFT_ACCELERATION
    }

    static inline Field projectForce(Field& P) {return P;}

    static inline void update_field(Field &P, Field &U, double ep)
    {
#ifndef USE_FFT_ACCELERATION
      double t0=usecond(); 
      U += P * ep;
      double t1=usecond();
      double total_time = (t1-t0)/1e6;
      std::cout << GridLogIntegrator << "Total time for updating field (s)       : " << total_time << std::endl; 
#else
      // FFT transform P(x) -> P(p)
      // divide by (M^2+p^2)  M external parameter (how to pass?)
      // P'(p) = P(p)/(M^2+p^2)
      // Transform back -> P'(x)
      // U += P'(x)*ep

      Field Pp(U._grid), P_FFT(U._grid);     
      static ComplexField p2(U._grid);
      RealD M = FFT_MASS;
      
      FFT theFFT((GridCartesian*)U._grid);
      theFFT.FFT_all_dim(Pp, P, FFT::forward);

      static bool first_call = true;
      if (first_call)
      {
        // avoid recomputing
        MomentumSpacePropagator(p2, M);
        first_call = false;
      }
      Pp *= p2;
      theFFT.FFT_all_dim(P_FFT, Pp, FFT::backward);
      U += P_FFT * ep;

#endif //USE_FFT_ACCELERATION
    }

    static inline RealD FieldSquareNorm(Field &U)
    {
#ifndef USE_FFT_ACCELERATION
      return (TensorRemove(sum(trace(U * U))).real());
#else
      // In case of Fourier acceleration we have to:
      // compute U(p)*U(p)/(M^2+p^2))   Parseval theorem
      // 1 FFT needed U(x) -> U(p)
      // M to be passed

      FFT theFFT((GridCartesian*)U._grid);
      Field Up(U._grid);

      theFFT.FFT_all_dim(Up, U, FFT::forward);
      RealD M = FFT_MASS;
      ComplexField p2(U._grid);
      MomentumSpacePropagator(p2, M);
      Field Up2 = Up * p2;
      // from the definition of the DFT we need to divide by the volume
      return (-TensorRemove(sum(trace(adj(Up) * Up2))).real() / U._grid->gSites());
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
