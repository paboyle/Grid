#ifndef  GRID_QCD_DHOP_H
#define  GRID_QCD_DHOP_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper routines that implement Wilson stencil for a single site.
    // Common to both the WilsonFermion and WilsonFermion5D
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Impl> class WilsonKernels : public FermionOperator<Impl> { 
    public:

     INHERIT_IMPL_TYPES(Impl);
     typedef FermionOperator<Impl> Base;
     
    public:
     void DiracOptDhopSite(StencilImpl &st,DoubledGaugeField &U,
			   std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			   int sF,int sU,const FermionField &in, FermionField &out);
      
     void DiracOptDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in,FermionField &out);

     void DiracOptDhopDir(StencilImpl &st,DoubledGaugeField &U,
			  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			  int sF,int sU,const FermionField &in, FermionField &out,int dirdisp,int gamma);

     void DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in, FermionField &out,uint64_t *);

     void DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
			       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			       int sF,int sU,const FermionField &in, FermionField &out);

     void DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
				  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
				  int sF,int sU,const FermionField &in, FermionField &out);

     WilsonKernels(const ImplParams &p= ImplParams());
     
    };

  }
}
#endif
