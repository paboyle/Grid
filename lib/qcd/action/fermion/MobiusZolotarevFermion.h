#ifndef  GRID_QCD_MOBIUS_ZOLOTAREV_FERMION_H
#define  GRID_QCD_MOBIUS_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class MobiusZolotarevFermion : public CayleyFermion5D
    {
    public:

      virtual void   Instantiatable(void) {};
      // Constructors
       MobiusZolotarevFermion(LatticeGaugeField &_Umu,
			      GridCartesian         &FiveDimGrid,
			      GridRedBlackCartesian &FiveDimRedBlackGrid,
			      GridCartesian         &FourDimGrid,
			      GridRedBlackCartesian &FourDimRedBlackGrid,
			      RealD _mass,RealD _M5,
			      RealD b, RealD c,
			      RealD lo, RealD hi) : 
      
      CayleyFermion5D(_Umu,
		      FiveDimGrid,
		      FiveDimRedBlackGrid,
		      FourDimGrid,
		      FourDimRedBlackGrid,_mass,_M5)

      {
	RealD eps = lo/hi;

	Approx::zolotarev_data *zdata = Approx::zolotarev(eps,this->Ls,0);
	assert(zdata->n==this->Ls);

	std::cout << "MobiusZolotarevFermion (b="<<b<<",c="<<c<<") with Ls= "<<Ls<<" Zolotarev range ["<<lo<<","<<hi<<"]"<<std::endl;
	
	// Call base setter
	this->CayleyFermion5D::SetCoefficientsZolotarev(hi,zdata,b,c);
 
	Approx::zolotarev_free(zdata);
      }

    };

  }
}

#endif
