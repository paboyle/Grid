#ifndef OVERLAP_WILSON_CONTFRAC_TANH_FERMION_H
#define OVERLAP_WILSON_CONTFRAC_TANH_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class OverlapWilsonContFracTanhFermion : public ContinuedFractionFermion5D
    {
    public:

      virtual void   Instantiatable(void){};
      // Constructors
    OverlapWilsonContFracTanhFermion(LatticeGaugeField &_Umu,
				     GridCartesian         &FiveDimGrid,
				     GridRedBlackCartesian &FiveDimRedBlackGrid,
				     GridCartesian         &FourDimGrid,
				     GridRedBlackCartesian &FourDimRedBlackGrid,
				     RealD _mass,RealD _M5,
				     RealD scale) :
      
      // b+c=scale, b-c = 0 <=> b =c = scale/2
      ContinuedFractionFermion5D(_Umu,
				 FiveDimGrid,
				 FiveDimRedBlackGrid,
				 FourDimGrid,
				 FourDimRedBlackGrid,_mass,_M5)
	{
	  assert((Ls&0x1)==1); // Odd Ls required
	  int nrational=Ls-1;// Even rational order
	  Approx::zolotarev_data *zdata = Approx::higham(1.0,nrational);// eps is ignored for higham
	  SetCoefficientsTanh(zdata,scale);
	  Approx::zolotarev_free(zdata);
	}
    };
  }
}
#endif
