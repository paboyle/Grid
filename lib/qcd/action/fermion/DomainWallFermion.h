#ifndef  GRID_QCD_DOMAIN_WALL_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    class DomainWallFermion : public CayleyFermion5D
    {
    public:

      virtual void   Instantiatable(void) {};
      // Constructors
      DomainWallFermion(LatticeGaugeField &_Umu,
			GridCartesian         &FiveDimGrid,
			GridRedBlackCartesian &FiveDimRedBlackGrid,
			GridCartesian         &FourDimGrid,
			GridRedBlackCartesian &FourDimRedBlackGrid,
			RealD _mass,RealD _M5) : 

      CayleyFermion5D(_Umu,
		      FiveDimGrid,
		      FiveDimRedBlackGrid,
		      FourDimGrid,
		      FourDimRedBlackGrid,_mass,_M5)

      {
	RealD eps = 1.0;

	Approx::zolotarev_data *zdata = Approx::grid_higham(eps,this->Ls);// eps is ignored for higham
	assert(zdata->n==this->Ls);
	
	std::cout << "DomainWallFermion with Ls="<<Ls<<std::endl;
	// Call base setter
	this->CayleyFermion5D::SetCoefficientsTanh(zdata,1.0,0.0);
 
      }

    };

  }
}

#endif
