#ifndef OVERLAP_WILSON_CONTFRAC_ZOLOTAREV_FERMION_H
#define OVERLAP_WILSON_CONTFRAC_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class OverlapWilsonContFracZolotarevFermion : public ContinuedFractionFermion5D<Impl>
    {
    public:
#include <qcd/action/fermion/FermionImplTypedefs.h>
    public:

      virtual void   Instantiatable(void){};
      // Constructors
    OverlapWilsonContFracZolotarevFermion(GaugeField &_Umu,
					  GridCartesian         &FiveDimGrid,
					  GridRedBlackCartesian &FiveDimRedBlackGrid,
					  GridCartesian         &FourDimGrid,
					  GridRedBlackCartesian &FourDimRedBlackGrid,
					  RealD _mass,RealD _M5,
					  RealD lo,RealD hi):
      
      // b+c=scale, b-c = 0 <=> b =c = scale/2
      ContinuedFractionFermion5D<Impl>(_Umu,
				 FiveDimGrid,
				 FiveDimRedBlackGrid,
				 FourDimGrid,
				 FourDimRedBlackGrid,_mass,_M5)
	{
	  assert((this->Ls&0x1)==1); // Odd Ls required

	  int nrational=this->Ls;// Odd rational order
	  RealD eps = lo/hi;

	  Approx::zolotarev_data *zdata = Approx::zolotarev(eps,nrational,0);
	  this->SetCoefficientsZolotarev(hi,zdata);
	  Approx::zolotarev_free(zdata);

	}
    };
  }
}
#endif
