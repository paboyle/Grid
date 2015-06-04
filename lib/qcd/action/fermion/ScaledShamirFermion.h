#ifndef  GRID_QCD_SCALED_SHAMIR_FERMION_H
#define  GRID_QCD_SCALED_SHAMIR_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class ScaledShamirFermion : public MobiusFermion
    {
    public:

      // Constructors
    ScaledShamirFermion(LatticeGaugeField &_Umu,
			GridCartesian         &FiveDimGrid,
			GridRedBlackCartesian &FiveDimRedBlackGrid,
			GridCartesian         &FourDimGrid,
			GridRedBlackCartesian &FourDimRedBlackGrid,
			RealD _mass,RealD _M5,
			RealD scale) :
      
      // b+c=scale, b-c = 1 <=> 2b = scale+1; 2c = scale-1
      MobiusFermion(_Umu,
		    FiveDimGrid,
		    FiveDimRedBlackGrid,
		    FourDimGrid,
		    FourDimRedBlackGrid,_mass,_M5,0.5*(scale+1.0),0.5*(scale-1.0))
      {
      }

    };

  }
}

#endif
