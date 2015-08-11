    public:
      typedef typename Impl::Simd                           Simd;
      typedef typename Impl::FermionField           FermionField;
      typedef typename Impl::GaugeLinkField       GaugeLinkField;
      typedef typename Impl::GaugeField               GaugeField;
      typedef typename Impl::DoubledGaugeField DoubledGaugeField;
      typedef typename Impl::SiteSpinor               SiteSpinor;
      typedef typename Impl::SiteHalfSpinor       SiteHalfSpinor;
      typedef typename Impl::Compressor               Compressor;
      typedef WilsonKernels<Impl> Kernels;
