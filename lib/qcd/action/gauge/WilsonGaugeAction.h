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
	return WilsonLoops<MatrixField,GaugeField>::sumPlaquette(U);
      };
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {
	//FIXME loop on directions
	MatrixField dSdU_mu(U._grid);
	WilsonLoops<MatrixField,GaugeField>::Staple(dSdU_mu,U,0);
      };
      virtual void  staple(const GaugeField &stap,GaugeField & U) {
	//FIXME loop on directions
	MatrixField stap_mu(U._grid);
	WilsonLoops<MatrixField,GaugeField>::Staple(stap_mu,U,0);
      };
    };
    
  }
}

#endif
