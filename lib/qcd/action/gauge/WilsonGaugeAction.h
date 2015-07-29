#ifndef QCD_WILSON_GAUGE_ACTION_H
#define QCD_WILSON_GAUGE_ACTION_H

namespace Grid{
  namespace QCD{
    
    ////////////////////////////////////////////////////////////////////////
    // Wilson Gauge Action .. should I template the Nc etc..
    ////////////////////////////////////////////////////////////////////////
    template<class GaugeField,class MatrixField>
      class WilsonGaugeAction : public Action<GaugeField> {
    private:
      RealD beta;
    public:
    WilsonGaugeAction(RealD b):beta(b){};
      
      virtual void init(const GaugeField &U, GridParallelRNG& pRNG) {};
      
      virtual RealD S(const GaugeField &U) {
	RealD plaq = WilsonLoops<MatrixField,GaugeField>::avgPlaquette(U);
	std::cout<<GridLogMessage << "Plaq : "<<plaq << "\n";
	RealD vol = U._grid->gSites();
	RealD action=beta*(1.0 -plaq)*(Nd*(Nd-1.0))*vol*0.5;
	std::cout << GridLogMessage << "WilsonGauge action "<<action<<std::endl;
	return action;
      };
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

	//not optimal implementation FIXME
	//extend Ta to include Lorentz indexes
	RealD factor = 0.5*beta/RealD(Nc);
	MatrixField dSdU_mu(U._grid);
	MatrixField Umu(U._grid);
	for (int mu=0; mu < Nd; mu++){
	  Umu = PeekIndex<LorentzIndex>(U,mu);
	  // Staple in direction mu
	  WilsonLoops<MatrixField,GaugeField>::Staple(dSdU_mu,U,mu);
	  dSdU_mu = Ta(Umu*adj(dSdU_mu))*factor;
	  pokeLorentz(dSdU, dSdU_mu, mu);
	}
      };
    };
    
  }
}

#endif
