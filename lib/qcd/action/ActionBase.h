#ifndef QCD_ACTION_BASE
#define QCD_ACTION_BASE
namespace Grid {
namespace QCD{

template<class GaugeField>
class Action { 

 public:
  // Boundary conditions? // Heatbath?
  virtual void  refresh(const GaugeField &U, GridParallelRNG& pRNG) = 0;// refresh pseudofermions
  virtual RealD S    (const GaugeField &U)                        = 0;  // evaluate the action
  virtual void  deriv(const GaugeField &U,GaugeField & dSdU )     = 0;  // evaluate the action derivative
  virtual ~Action() {};
};

// Could derive PseudoFermion action with a PF field, FermionField, and a Grid; implement refresh
/*
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
*/

template<class GaugeField> struct ActionLevel{
public:
   
  typedef Action<GaugeField>*  ActPtr; // now force the same colours as the rest of the code

  int multiplier;

  std::vector<ActPtr> actions;

  ActionLevel(int mul = 1) : multiplier(mul) {
    assert (mul > 0);
  };
   
  void push_back(ActPtr ptr){
    actions.push_back(ptr);
  }
};

template<class GaugeField> using ActionSet = std::vector<ActionLevel< GaugeField > >;


}}
#endif
