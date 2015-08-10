#ifndef  OVERLAP_WILSON_CAYLEY_ZOLOTAREV_FERMION_H
#define  OVERLAP_WILSON_CAYLEY_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class OverlapWilsonCayleyZolotarevFermion : public MobiusZolotarevFermion<Impl>
    {
    public:
#include <qcd/action/fermion/FermionImplTypedefs.h>
    public:

      // Constructors

    OverlapWilsonCayleyZolotarevFermion(GaugeField &_Umu,
					GridCartesian         &FiveDimGrid,
					GridRedBlackCartesian &FiveDimRedBlackGrid,
					GridCartesian         &FourDimGrid,
					GridRedBlackCartesian &FourDimRedBlackGrid,
					RealD _mass,RealD _M5,
					RealD lo, RealD hi) : 
      // b+c=1.0, b-c = 0 <=> b =c = 1/2
      MobiusZolotarevFermion<Impl>(_Umu,
			     FiveDimGrid,
			     FiveDimRedBlackGrid,
			     FourDimGrid,
			     FourDimRedBlackGrid,_mass,_M5,0.5,0.5,lo,hi)

      {}

    };

  }
}

#endif
