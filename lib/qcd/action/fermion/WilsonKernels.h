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
#if 0
     void DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in, FermionField &out,uint64_t *p){
       DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
     }
#endif
// doesn't seem to work with Gparity at the moment
#undef HANDOPT
#if 1
     void DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
			       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			       int sF,int sU,const FermionField &in, FermionField &out);

     void DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
				  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
				  int sF,int sU,const FermionField &in, FermionField &out);
#endif

     WilsonKernels(const ImplParams &p= ImplParams());
     
    };

  }
}
#endif
