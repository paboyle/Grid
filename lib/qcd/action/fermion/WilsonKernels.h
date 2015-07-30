#ifndef  GRID_QCD_DHOP_H
#define  GRID_QCD_DHOP_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper classes that implement Wilson stencil for a single site.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Generic version works for any Nc and with extra flavour indices
    //    namespace DiracOpt {

      // These ones will need to be package intelligently. WilsonType base class
      // for use by DWF etc..
      void DiracOptDhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
			    std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			    int sF,int sU,const LatticeFermion &in, LatticeFermion &out);
      void DiracOptDhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
			       std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			       int sF,int sU,const LatticeFermion &in, LatticeFermion &out);
      void DiracOptDhopDir(CartesianStencil &st,LatticeDoubledGaugeField &U,
			   std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
			   int sF,int sU,const LatticeFermion &in, LatticeFermion &out,int dirdisp,int gamma);
      
      //  };

      // Hand unrolled for Nc=3, one flavour
      //    namespace DiracOptHand {
      // These ones will need to be package intelligently. WilsonType base class
      // for use by DWF etc..

      void DiracOptHandDhopSite(CartesianStencil &st,LatticeDoubledGaugeField &U,
				std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
				int sF,int sU,const LatticeFermion &in, LatticeFermion &out);
      void DiracOptHandDhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
				   std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
				   int sF,int sU,const LatticeFermion &in, LatticeFermion &out);

      //    };
  

    void DiracOptHandDhopSiteDag(CartesianStencil &st,LatticeDoubledGaugeField &U,
				 std::vector<vHalfSpinColourVector,alignedAllocator<vHalfSpinColourVector> >  &buf,
				 int sF,int sU,const LatticeFermion &in, LatticeFermion &out);

  }
}
#endif
