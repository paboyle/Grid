#ifndef  GRID_QCD_DOMAIN_WALL_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class ScaledShamirFermion : public CayleyFermion5D
    {
    public:

      // Constructors
      ScaledShamirFermion(LatticeGaugeField &_Umu,
			  GridCartesian         &FiveDimGrid,
			  GridRedBlackCartesian &FiveDimRedBlackGrid,
			  GridCartesian         &FourDimGrid,
			  GridRedBlackCartesian &FourDimRedBlackGrid,
			  RealD _mass,RealD _M5, RealD scale) : 
      
      CayleyFermion5D(_Umu,
		      FiveDimGrid,
		      FiveDimRedBlackGrid,
		      FourDimGrid,
		      FourDimRedBlackGrid,_mass,_M5,
		      RealD b, 
		      RealD c)

      {
	RealD eps = 1.0;

	Approx::zolotarev_data *zdata = Approx::grid_higham(eps,this->Ls);// eps is ignored for higham
	assert(zdata->n==this->Ls);
	
	//b+c = scale;
	//b-c = 1
	//b   = 0.5(scale+1);
	//c   = 0.5(scale-1);
	
	// Call base setter
	this->CayleyFermion5D::SetCoefficients(1.0,zdata,0.5*(scale+1.0),0.5*(scale-1.0));

       }

    };

  }
}

#endif
