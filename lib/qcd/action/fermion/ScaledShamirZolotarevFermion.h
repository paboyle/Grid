#ifndef  GRID_QCD_SCALED_SHAMIR_ZOLOTAREV_FERMION_H
#define  GRID_QCD_SCALED_SHAMIR_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class ScaledShamirZolotarevFermion : public MobiusZolotarevFermion
    {
    public:

      // Constructors


    ScaledShamirZolotarevFermion(LatticeGaugeField &_Umu,
				 GridCartesian         &FiveDimGrid,
				 GridRedBlackCartesian &FiveDimRedBlackGrid,
				 GridCartesian         &FourDimGrid,
				 GridRedBlackCartesian &FourDimRedBlackGrid,
				 RealD _mass,RealD _M5,
				 RealD scale,
				 RealD lo, RealD hi) : 
      
      MobiusZolotarevFermion(_Umu,
			       FiveDimGrid,
			       FiveDimRedBlackGrid,
			       FourDimGrid,
			       FourDimRedBlackGrid,_mass,_M5,0.5*(scale+1.0),0.5*(scale-1.0),lo,hi)

      {}

    };

  }
}

#endif
