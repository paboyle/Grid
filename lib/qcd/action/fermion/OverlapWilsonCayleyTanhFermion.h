#ifndef OVERLAP_WILSON_CAYLEY_TANH_FERMION_H
#define OVERLAP_WILSON_CAYLEY_TANH_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class OverlapWilsonCayleyTanhFermion : public MobiusFermion
    {
    public:

      // Constructors
    OverlapWilsonCayleyTanhFermion(LatticeGaugeField &_Umu,
				   GridCartesian         &FiveDimGrid,
				   GridRedBlackCartesian &FiveDimRedBlackGrid,
				   GridCartesian         &FourDimGrid,
				   GridRedBlackCartesian &FourDimRedBlackGrid,
				   RealD _mass,RealD _M5,
				   RealD scale) :
      
      // b+c=scale, b-c = 0 <=> b =c = scale/2
      MobiusFermion(_Umu,
		    FiveDimGrid,
		    FiveDimRedBlackGrid,
		    FourDimGrid,
		    FourDimRedBlackGrid,_mass,_M5,0.5*scale,0.5*scale)
	{
	}
    };
  }
}
#endif
