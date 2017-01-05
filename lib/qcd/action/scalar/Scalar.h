#ifndef QCD_SCALAR_ACTION_H
#define QCD_SCALAR_ACTION_H

#define INHERIT_SIMPL_TYPES(Impl)\
typedef typename Impl::SiteScalar      SiteScalar;	 \
typedef typename Impl::SiteSpinor      SiteSpinor;	 \
typedef typename Impl::SitePropagator  SitePropagator;   \
typedef typename Impl::ScalarField     ScalarField;	 \
typedef typename Impl::FermionField    FermionField;	 \
typedef typename Impl::PropagatorField PropagatorField;  \
typedef typename Impl::StencilImpl     StencilImpl;

namespace Grid{
namespace QCD{
  // Scalar implementation class ///////////////////////////////////////////////
  // FIXME: it is not very nice to have the FImpl aliases
  template <class S,
            class Representation = FundamentalRep<1>,
            class _Coeff_t = RealD>
  class ScalarImpl:
    public PeriodicGaugeImpl<GaugeImplTypes<S, Representation::Dimension>>
  {
  public:
    static constexpr unsigned int rDim = Representation::Dimension;
  public:
    // gauge types
    typedef PeriodicGaugeImpl<GaugeImplTypes<S, rDim>> Gimpl;
    INHERIT_GIMPL_TYPES(Gimpl);
    // site types
    // (using classes instead of aliases to allow for partial specialisation)
    template <typename vtype, unsigned int d>
    class iImplScalar
    {
    public:
      typedef iScalar<iScalar<iVector<vtype, d>>> type;
    };
    template <typename vtype>
    class iImplScalar<vtype, 1>
    {
    public:
      typedef iScalar<iScalar<iScalar<vtype>>> type;
    };
    template <typename vtype, unsigned int d>
    class iImplPropagator
    {
    public:
      typedef iScalar<iScalar<iMatrix<vtype, d>>> type;
    };
    template <typename vtype>
    class iImplPropagator<vtype, 1>
    {
    public:
      typedef iScalar<iScalar<iScalar<vtype>>> type;
    };
    // type aliases
    typedef typename iImplScalar<S, rDim>::type      SiteScalar;
    typedef SiteScalar                               SiteSpinor;
    typedef typename iImplPropagator<S, rDim>::type  SitePropagator;
    typedef Lattice<SiteScalar>                      ScalarField;
    typedef ScalarField                              FermionField;
    typedef Lattice<SitePropagator>                  PropagatorField;
    typedef CartesianStencil<SiteScalar, SiteScalar> StencilImpl;
  };
  
  // single scalar implementation
  typedef ScalarImpl<vComplex> ScalarImplR;
  
  // Scalar action /////////////////////////////////////////////////////////////
  template <typename SImpl>
  class Scalar:
    public CheckerBoardedSparseMatrixBase<typename SImpl::ScalarField>,
    public SImpl
  {
  public:
    INHERIT_GIMPL_TYPES(SImpl);
    INHERIT_SIMPL_TYPES(SImpl);
  public:
    // constructor
    Scalar(GaugeField &_Umu, GridCartesian &Sgrid, GridRedBlackCartesian &Hgrid,
           RealD _mass)
    : _grid(&Sgrid)
    , _cbgrid(&Hgrid)
    , mass(_mass)
    , Lebesgue(_grid)
    , LebesgueEvenOdd(_cbgrid)
    , Umu(&Sgrid)
    , UmuEven(&Hgrid)
    , UmuOdd(&Hgrid)
    {
      Umu = _Umu;
      pickCheckerboard(Even, UmuEven, Umu);
      pickCheckerboard(Odd, UmuOdd, Umu);
    }
    // grid access
    virtual GridBase *RedBlackGrid(void) {return _grid;}
    // half checkerboard operations
    // FIXME: do implementation
    virtual void Meooe(const ScalarField &in, ScalarField &out)
    {
      assert(0);
    }
    virtual void Mooee(const ScalarField &in, ScalarField &out)
    {
      assert(0);
    }
    virtual void MooeeInv(const ScalarField &in, ScalarField &out)
    {
      assert(0);
    }
    virtual void MeooeDag(const ScalarField &in, ScalarField &out)
    {
      assert(0);
    }
    virtual void MooeeDag(const ScalarField &in, ScalarField &out)
    {
      assert(0);
    }
    virtual void MooeeInvDag(const ScalarField &in, ScalarField &out)
    {
      assert(0);
    }
    // free propagators
    static void MomentumSpacePropagator(ScalarField &out, RealD m);
    static void FreePropagator(const ScalarField &in, ScalarField &out,
                               const ScalarField &momKernel);
    static void FreePropagator(const ScalarField &in, ScalarField &out, RealD m);
  public:
    RealD mass;
    
    GridBase *_grid;
    GridBase *_cbgrid;
    
    // Defines the stencils for even and odd
    StencilImpl Stencil;
    StencilImpl StencilEven;
    StencilImpl StencilOdd;
    
    // Copy of the gauge field, with even and odd subsets
    GaugeField Umu;
    GaugeField UmuEven;
    GaugeField UmuOdd;
    
    LebesgueOrder Lebesgue;
    LebesgueOrder LebesgueEvenOdd;
  };
  
  template <typename SImpl>
  void Scalar<SImpl>::MomentumSpacePropagator(ScalarField &out, RealD m)
  {
    GridBase           *grid = out._grid;
    ScalarField        kmu(grid), one(grid);
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
  
  template <typename SImpl>
  void Scalar<SImpl>::FreePropagator(const ScalarField &in, ScalarField &out,
                                     const ScalarField &FTKernel)
  {
    FFT         fft((GridCartesian *)in._grid);
    ScalarField inFT(in._grid);
    
    fft.FFT_all_dim(inFT, in, FFT::forward);
    inFT = inFT*FTKernel;
    fft.FFT_all_dim(out, inFT, FFT::backward);
  }
  
  template <typename SImpl>
  void Scalar<SImpl>::FreePropagator(const ScalarField &in, ScalarField &out,
                                     RealD m)
  {
    ScalarField FTKernel(in._grid);
    
    MomentumSpacePropagator(FTKernel, m);
    FreePropagator(in, out, FTKernel);
  }
  
  template <class SImpl>
  void ScalarToProp(typename SImpl::PropagatorField &p,
                    const typename SImpl::ScalarField &s,
                    const int c)
  {
    for(int i = 0; i < SImpl::rDim; ++i)
    {
      pokeColour(p, peekColour(s, i), i);
    }
  }
  
  template <class SImpl>
  void PropToScalar(typename SImpl::ScalarField &s,
                    const typename SImpl::PropagatorField &p,
                    const int c)
  {
    for(int i = 0; i < SImpl::rDim; ++i)
    {
      pokeColour(s, peekColour(p, i), i);
    }
  }
}}

#endif
