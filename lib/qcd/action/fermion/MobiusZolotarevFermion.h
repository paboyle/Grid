#ifndef  GRID_QCD_MOBIUS_ZOLOTAREV_FERMION_H
#define  GRID_QCD_MOBIUS_ZOLOTAREV_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class MobiusZolotarevFermion : public CayleyFermion5D
    {
    public:

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

	Approx::zolotarev_data *zdata = Approx::grid_zolotarev(eps,this->Ls,0);// eps is ignored for higham
	assert(zdata->n==this->Ls);

	std::cout << "MobiusZolotarevFermion (b="<<b<<",c="<<c<<") with Ls= "<<Ls<<" Zolotarev range ["<<lo<<","<<hi<<"]"<<std::endl;
	
	// Call base setter
	this->CayleyFermion5D::SetCoefficients(1.0,zdata,b,c);
 
      }

    };

  }
}

#endif
