#ifndef OVERLAP_WILSON_CAYLEY_TANH_FERMION_H
#define OVERLAP_WILSON_CAYLEY_TANH_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class OverlapWilsonContFracZolotarevFermion : public ContinuedFractionFermion5D
    {
    public:

      virtual void   Instantiatable(void){};
      // Constructors
    OverlapWilsonContFracZolotarevFermion(LatticeGaugeField &_Umu,
					  GridCartesian         &FiveDimGrid,
					  GridRedBlackCartesian &FiveDimRedBlackGrid,
					  GridCartesian         &FourDimGrid,
					  GridRedBlackCartesian &FourDimRedBlackGrid,
					  RealD _mass,RealD _M5,
					  RealD lo,RealD hi):
      
      // b+c=scale, b-c = 0 <=> b =c = scale/2
      ContinuedFractionFermion5D(_Umu,
				 FiveDimGrid,
				 FiveDimRedBlackGrid,
				 FourDimGrid,
				 FourDimRedBlackGrid,_mass)
	{
	  assert((Ls&0x1)==1); // Odd Ls required

	  int nrational=Ls-1;// Even rational order
	  RealD eps = lo/hi;

	  Approx::zolotarev_data *zdata = Approx::grid_zolotarev(eps,nrational,0);

	  SetCoefficientsZolotarev(hi,zdata);

	}
    };
  }
}
#endif
