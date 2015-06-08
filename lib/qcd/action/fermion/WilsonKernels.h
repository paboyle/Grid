#ifndef  GRID_QCD_DHOP_H
#define  GRID_QCD_DHOP_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper classes that implement Wilson stencil for a single site.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Generic version works for any Nc and with extra flavour indices
    class DiracOpt {
    public:
      // These ones will need to be package intelligently. WilsonType base class
      // for use by DWF etc..
      static void DhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
			   std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			   int sF,int sU,const LatticeFermion &in, LatticeFermion &out);
      static void DhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
			      std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			      int sF,int sU,const LatticeFermion &in, LatticeFermion &out);
      static void DhopDir(CartesianStencil &st,LatticeDoubledGaugeField &U,
			   std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			  int sF,int sU,const LatticeFermion &in, LatticeFermion &out,int dirdisp);

    };

    // Hand unrolled for Nc=3, one flavour
    class DiracOptHand {
    public:
      // These ones will need to be package intelligently. WilsonType base class
      // for use by DWF etc..
      static void DhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
			   std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			   int sF,int sU,const LatticeFermion &in, LatticeFermion &out);
      static void DhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
			      std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			      int sF,int sU,const LatticeFermion &in, LatticeFermion &out);

    };

  }
}
#endif
