#ifndef  GRID_QCD_DOMAIN_WALL_FERMION_H
#define  GRID_QCD_DOMAIN_WALL_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class DomainWallFermion : public CayleyFermion5D<Impl>
    {
    public:
#include <qcd/action/fermion/FermionImplTypedefs.h>
    public:

      virtual void   Instantiatable(void) {};
      // Constructors
      DomainWallFermion(GaugeField &_Umu,
			GridCartesian         &FiveDimGrid,
			GridRedBlackCartesian &FiveDimRedBlackGrid,
			GridCartesian         &FourDimGrid,
			GridRedBlackCartesian &FourDimRedBlackGrid,
			RealD _mass,RealD _M5) : 

      CayleyFermion5D<Impl>(_Umu,
			    FiveDimGrid,
			    FiveDimRedBlackGrid,
			    FourDimGrid,
			    FourDimRedBlackGrid,_mass,_M5)

      {
	RealD eps = 1.0;

	Approx::zolotarev_data *zdata = Approx::higham(eps,this->Ls);// eps is ignored for higham
	assert(zdata->n==this->Ls);
	
	std::cout<<GridLogMessage << "DomainWallFermion with Ls="<<this->Ls<<std::endl;
	// Call base setter
	this->SetCoefficientsTanh(zdata,1.0,0.0);

	Approx::zolotarev_free(zdata);
      }

    };

  }
}

#endif
