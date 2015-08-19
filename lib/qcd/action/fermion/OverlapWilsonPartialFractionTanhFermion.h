#ifndef OVERLAP_WILSON_PARTFRAC_TANH_FERMION_H
#define OVERLAP_WILSON_PARTFRAC_TANH_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class OverlapWilsonPartialFractionTanhFermion : public PartialFractionFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      virtual void   Instantiatable(void){};
      // Constructors
    OverlapWilsonPartialFractionTanhFermion(GaugeField &_Umu,
					    GridCartesian         &FiveDimGrid,
					    GridRedBlackCartesian &FiveDimRedBlackGrid,
					    GridCartesian         &FourDimGrid,
					    GridRedBlackCartesian &FourDimRedBlackGrid,
					    RealD _mass,RealD _M5,
					    RealD scale,const ImplParams &p= ImplParams()) :
      
      // b+c=scale, b-c = 0 <=> b =c = scale/2
      PartialFractionFermion5D<Impl>(_Umu,
				     FiveDimGrid,
				     FiveDimRedBlackGrid,
				     FourDimGrid,
				     FourDimRedBlackGrid,_mass,_M5,p)
	{
	  assert((this->Ls&0x1)==1); // Odd Ls required
	  int nrational=this->Ls-1;// Even rational order
	  Approx::zolotarev_data *zdata = Approx::higham(1.0,nrational);// eps is ignored for higham
	  this->SetCoefficientsTanh(zdata,scale);
	  Approx::zolotarev_free(zdata);
	}
    };
  }
}
#endif
