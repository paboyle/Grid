#ifndef QCD_GAUGE_ACTION_BASE
#define QCD_GAUGE_ACTION_BASE
namespace Grid {
namespace QCD{

template<class GaugeField>
class GaugeActionBase { // derive this from TermInAction?

 public:
  virtual RealD S(const GaugeField &U)                        = 0;  // evaluate the action
  virtual void  deriv(const GaugeField &U,GaugeField & dSdU ) = 0;  // evaluate the action derivative
  virtual void  staple(const GaugeField &U,GaugeField & dSdU ) = 0;  // evaluate the action derivative
  virtual void refresh(const GaugeField & ) {}; 
  // Boundary conditions?
  // Heatbath?
  virtual ~GaugeActionBase() {};
};

}}
#endif
