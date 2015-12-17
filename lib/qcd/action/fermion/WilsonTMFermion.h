#ifndef  GRID_QCD_WILSON_TM_FERMION_H
#define  GRID_QCD_WILSON_TM_FERMION_H

#include <Grid.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class WilsonTMFermion : public WilsonFermion<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:
     RealD mu; // TwistedMass parameter

      virtual void   Instantiatable(void) {};
      // Constructors
      WilsonTMFermion(GaugeField &_Umu,
		    GridCartesian         &Fgrid,
		    GridRedBlackCartesian &Hgrid, 
		    RealD _mass,
		    RealD _mu,
		    const ImplParams &p= ImplParams()
		      ) :
	WilsonFermion<Impl>(_Umu,
			    Fgrid,
			    Hgrid,
			    _mass,p)

      {
	mu = _mu;
      }

    };
    // allow override for twisted mass and clover
    virtual void Mooee(const FermionField &in, FermionField &out) ;
    virtual void MooeeDag(const FermionField &in, FermionField &out) ;
    virtual void MooeeInv(const FermionField &in, FermionField &out) ;
    virtual void MooeeInvDag(const FermionField &in, FermionField &out) ;
    

  }
}

#endif
