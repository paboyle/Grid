#ifndef QCD_WILSON_GAUGE_ACTION_H
#define QCD_WILSON_GAUGE_ACTION_H

////////////////////////////////////////////////////////////////////////
// Wilson Gauge Action .. should I template the Nc etc..
////////////////////////////////////////////////////////////////////////
template<class GaugeField,class MatrixField>
class WilsonGaugeAction : public GaugeActionBase<GaugeField> {
 public:
  virtual RealD S(const GaugeField &U) {
    return WilsonLoops<MatrixField,GaugeField>::sumPlaquette(U);
  };
  virtual RealD deriv(GaugeField &U,GaugeField & dSdU ) {
    WilsonLoops<MatrixField,GaugeField>::Staple(dSdU,U,mu);
  };
  virtual void  staple(const MatrixField &stap,GaugeField & U,int mu ) {
    WilsonLoops<MatrixField,GaugeField>::Staple(stap,U,mu);
  };
};


#endif

#endif
