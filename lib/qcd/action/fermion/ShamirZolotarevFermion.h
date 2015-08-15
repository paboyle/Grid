#ifndef  GRID_QCD_SHAMIR_ZOLOTAREV_FERMION_H
#define  GRID_QCD_SHAMIR_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class ShamirZolotarevFermion : public MobiusZolotarevFermion<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      // Constructors


    ShamirZolotarevFermion(GaugeField &_Umu,
			   GridCartesian         &FiveDimGrid,
			   GridRedBlackCartesian &FiveDimRedBlackGrid,
			   GridCartesian         &FourDimGrid,
			   GridRedBlackCartesian &FourDimRedBlackGrid,
			   RealD _mass,RealD _M5,
			   RealD lo, RealD hi,const ImplParams &p= ImplParams()) : 
      
      // b+c = 1; b-c = 1 => b=1, c=0
      MobiusZolotarevFermion<Impl>(_Umu,
				   FiveDimGrid,
				   FiveDimRedBlackGrid,
				   FourDimGrid,
				   FourDimRedBlackGrid,_mass,_M5,1.0,0.0,lo,hi,p)
      
      {}

    };

  }
}

#endif
