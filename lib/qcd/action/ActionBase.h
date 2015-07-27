#ifndef QCD_ACTION_BASE
#define QCD_ACTION_BASE
namespace Grid {
namespace QCD{

template<class GaugeField>
class Action { 

 public:
  virtual void  init (const GaugeField &U, GridParallelRNG& pRNG) = 0;  // 
  virtual RealD S    (const GaugeField &U)                        = 0;  // evaluate the action
  virtual void  deriv(const GaugeField &U,GaugeField & dSdU )     = 0;  // evaluate the action derivative
  virtual void  refresh(const GaugeField & ) {};                        // Default to no-op for actions with no internal fields
  // Boundary conditions?
  // Heatbath?
  virtual ~Action() {};
};

// Could derive PseudoFermion action with a PF field, FermionField, and a Grid; implement refresh
template<class GaugeField, class FermionField>
class PseudoFermionAction : public Action<GaugeField> {
 public:
  FermionField Phi;
  GridParallelRNG &pRNG;
  GridBase &Grid;

  PseudoFermionAction(GridBase &_Grid,GridParallelRNG &_pRNG) : Grid(_Grid), Phi(&_Grid), pRNG(_pRNG) {
  };

  virtual void refresh(const GaugeField &gauge) {
    gaussian(Phi,pRNG);
  };

};
}}
#endif
