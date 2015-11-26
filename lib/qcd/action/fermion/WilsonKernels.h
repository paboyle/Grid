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
#if defined(AVX512) || defined(IMCI)
     void DiracOptAsmDhopSite(CartesianStencil &st,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in, FermionField &out,uint64_t *);
#else
     //void DiracOptAsmDhopSite(CartesianStencil &st,DoubledGaugeField &U,
     void DiracOptAsmDhopSite(StencilImpl &st,DoubledGaugeField &U,
			      std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			      int sF,int sU,const FermionField &in, FermionField &out,uint64_t *p){
       DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
     }
#endif
#define HANDOPT
#ifdef HANDOPT
     void DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
			       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			       int sF,int sU,const FermionField &in, FermionField &out);

     void DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
				  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
				  int sF,int sU,const FermionField &in, FermionField &out);
#else

     void DiracOptHandDhopSite(StencilImpl &st,DoubledGaugeField &U,
			       std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			       int sF,int sU,const FermionField &in, FermionField &out)
     {
       DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
     }

     void DiracOptHandDhopSiteDag(StencilImpl &st,DoubledGaugeField &U,
				  std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
				  int sF,int sU,const FermionField &in, FermionField &out)
     {
       DiracOptDhopSiteDag(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
     }
#endif

     WilsonKernels(const ImplParams &p= ImplParams()) : Base(p) {};

    };

  }
}
#endif
