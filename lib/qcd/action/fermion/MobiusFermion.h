#ifndef  GRID_QCD_MOBIUS_FERMION_H
#define  GRID_QCD_MOBIUS_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class MobiusFermion : public CayleyFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      virtual void   Instantiatable(void) {};
      // Constructors
      MobiusFermion(GaugeField &_Umu,
		    GridCartesian         &FiveDimGrid,
		    GridRedBlackCartesian &FiveDimRedBlackGrid,
		    GridCartesian         &FourDimGrid,
		    GridRedBlackCartesian &FourDimRedBlackGrid,
		    RealD _mass,RealD _M5,
		    RealD b, RealD c,const ImplParams &p= ImplParams()) : 
      
      CayleyFermion5D<Impl>(_Umu,
			    FiveDimGrid,
			    FiveDimRedBlackGrid,
			    FourDimGrid,
			    FourDimRedBlackGrid,_mass,_M5,p)

      {
	RealD eps = 1.0;

	std::cout<<GridLogMessage << "MobiusFermion (b="<<b<<",c="<<c<<") with Ls= "<<this->Ls<<" Tanh approx"<<std::endl;
	Approx::zolotarev_data *zdata = Approx::higham(eps,this->Ls);// eps is ignored for higham
	assert(zdata->n==this->Ls);
	
	// Call base setter
	this->SetCoefficientsTanh(zdata,b,c);

	Approx::zolotarev_free(zdata);
 
      }

    };

  }
}

#endif
