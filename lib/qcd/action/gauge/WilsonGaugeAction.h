#ifndef QCD_WILSON_GAUGE_ACTION_H
#define QCD_WILSON_GAUGE_ACTION_H

namespace Grid{
  namespace QCD{
    
    ////////////////////////////////////////////////////////////////////////
    // Wilson Gauge Action .. should I template the Nc etc..
    ////////////////////////////////////////////////////////////////////////
    template<class GaugeField>
    class WilsonGaugeAction : public Action<GaugeField> {
    public:

      typedef LorentzScalar<GaugeField> GaugeLinkField;

    private:
      RealD beta;
    public:
    WilsonGaugeAction(RealD b):beta(b){};
      
      virtual void refresh(const GaugeField &U, GridParallelRNG& pRNG) {}; // noop as no pseudoferms
      
      virtual RealD S(const GaugeField &U) {
	RealD plaq = WilsonLoops<GaugeField>::avgPlaquette(U);
	RealD vol = U._grid->gSites();
	RealD action=beta*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5;
	return action;
      };

      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
	//not optimal implementation FIXME
	//extend Ta to include Lorentz indexes

	//RealD factor = 0.5*beta/RealD(Nc);
	RealD factor = 0.5*beta/RealD(Nc);

	GaugeLinkField Umu(U._grid);
	GaugeLinkField dSdU_mu(U._grid);
	for (int mu=0; mu < Nd; mu++){

	  Umu = PeekIndex<LorentzIndex>(U,mu);

	  // Staple in direction mu
	  WilsonLoops<GaugeField>::Staple(dSdU_mu,U,mu);
	  dSdU_mu = Ta(Umu*dSdU_mu)*factor;
	  PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
	}
      };
    };
    
  }
}

#endif
