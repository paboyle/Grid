#ifndef  GRID_QCD_MOBIUS_ZOLOTAREV_FERMION_H
#define  GRID_QCD_MOBIUS_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class MobiusZolotarevFermion : public CayleyFermion5D<Impl>
    {
    public:
#include <qcd/action/fermion/FermionImplTypedefs.h>
    public:

      virtual void   Instantiatable(void) {};
      // Constructors
       MobiusZolotarevFermion(GaugeField &_Umu,
			      GridCartesian         &FiveDimGrid,
			      GridRedBlackCartesian &FiveDimRedBlackGrid,
			      GridCartesian         &FourDimGrid,
			      GridRedBlackCartesian &FourDimRedBlackGrid,
			      RealD _mass,RealD _M5,
			      RealD b, RealD c,
			      RealD lo, RealD hi) : 
      
      CayleyFermion5D<Impl>(_Umu,
		      FiveDimGrid,
		      FiveDimRedBlackGrid,
		      FourDimGrid,
		      FourDimRedBlackGrid,_mass,_M5)

      {
	RealD eps = lo/hi;

	Approx::zolotarev_data *zdata = Approx::zolotarev(eps,this->Ls,0);
	assert(zdata->n==this->Ls);

	std::cout<<GridLogMessage << "MobiusZolotarevFermion (b="<<b<<",c="<<c<<") with Ls= "<<this->Ls<<" Zolotarev range ["<<lo<<","<<hi<<"]"<<std::endl;
	
	// Call base setter
	this->SetCoefficientsZolotarev(hi,zdata,b,c);
 
	Approx::zolotarev_free(zdata);
      }

    };

  }
}

#endif
