#ifndef QCD_ACTION_BASE
#define QCD_ACTION_BASE
namespace Grid {
namespace QCD{

template<class GaugeField>
class Action { 

 public:
  virtual void  init(const GaugeField &U, GridParallelRNG& pRNG) = 0;
  virtual RealD S(const GaugeField &U)                           = 0;  // evaluate the action
  virtual void  deriv(const GaugeField &U,GaugeField & dSdU )    = 0;  // evaluate the action derivative
  //virtual void  refresh(const GaugeField & ) {}                ; 
  // Boundary conditions?
  // Heatbath?
  virtual ~Action() {};
};

}}
#endif
