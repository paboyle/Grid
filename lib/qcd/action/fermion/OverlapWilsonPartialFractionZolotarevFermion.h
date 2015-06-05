#ifndef OVERLAP_WILSON_PARTFRAC_ZOLOTAREV_FERMION_H
#define OVERLAP_WILSON_PARTFRAC_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class OverlapWilsonPartialFractionZolotarevFermion : public PartialFractionFermion5D
    {
    public:

      virtual void   Instantiatable(void){};
      // Constructors
    OverlapWilsonPartialFractionZolotarevFermion(LatticeGaugeField &_Umu,
					  GridCartesian         &FiveDimGrid,
					  GridRedBlackCartesian &FiveDimRedBlackGrid,
					  GridCartesian         &FourDimGrid,
					  GridRedBlackCartesian &FourDimRedBlackGrid,
					  RealD _mass,RealD _M5,
					  RealD lo,RealD hi):
      
      // b+c=scale, b-c = 0 <=> b =c = scale/2
      PartialFractionFermion5D(_Umu,
			       FiveDimGrid,
			       FiveDimRedBlackGrid,
			       FourDimGrid,
			       FourDimRedBlackGrid,_mass,_M5)
	{
	  assert((Ls&0x1)==1); // Odd Ls required

	  int nrational=Ls;// Odd rational order
	  RealD eps = lo/hi;

	  Approx::zolotarev_data *zdata = Approx::zolotarev(eps,nrational,0);
	  SetCoefficientsZolotarev(hi,zdata);
	  Approx::zolotarev_free(zdata);

	}
    };
  }
}
#endif
