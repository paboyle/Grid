#ifndef  GRID_QCD_DHOP_H
#define  GRID_QCD_DHOP_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper routines that implement Wilson stencil for a single site.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Impl> 
    class WilsonKernels { 
    public:

#include <qcd/action/fermion/FermionImplTypedefs.h>

    public:
      static void 
	DiracOptDhopSite(CartesianStencil &st,DoubledGaugeField &U,
			 std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			 int sF,int sU,const FermionField &in, FermionField &out);
      
      static void 
	DiracOptDhopSiteDag(CartesianStencil &st,DoubledGaugeField &U,
			    std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			    int sF,int sU,const FermionField &in,FermionField &out);

      static void 
	DiracOptDhopDir(CartesianStencil &st,DoubledGaugeField &U,
			std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			int sF,int sU,const FermionField &in, FermionField &out,int dirdisp,int gamma);

      static void 
	DiracOptHandDhopSite(CartesianStencil &st,DoubledGaugeField &U,
			     std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
			     int sF,int sU,const FermionField &in, FermionField &out){
	DiracOptDhopSite(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
      }

      static void 
	DiracOptHandDhopSiteDag(CartesianStencil &st,DoubledGaugeField &U,
				std::vector<SiteHalfSpinor,alignedAllocator<SiteHalfSpinor> >  &buf,
				int sF,int sU,const FermionField &in, FermionField &out){
	DiracOptDhopSiteDag(st,U,buf,sF,sU,in,out); // will template override for Wilson Nc=3
      }

    };

  }
}
#endif
