//--------------------------------------------------------------------
/*! @file Integrator_base.h
 * @brief Declaration of classes for the abstract Molecular Dynamics integrator
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

class Action;
class RandNum;
class Observer;

typedef std::vector<Action*> ActionLevel;
typedef std::vector<ActionLevel> ActionSet;
typedef std::vector<Observer*> ObserverList;

/*! @brief Abstract base class for Molecular Dynamics management */

namespace Grid{
  namespace QCD{
    
    class Integrator{
    private:
      virtual void update_P(int lv,double ep) = 0;
      virtual void update_U(double ep) = 0;
      
      virtual void register_observers() = 0;
      virtual void notify_observers() = 0;
      
    public:
      virtual ~Integrator(){}
      virtual void init(const LatticeColourMatrix&,
			const GridParallelRNG& RNG)=0;
      virtual double S()const =0;
      virtual void integrate(int level) =0;
      virtual const LatticeColourMatrix get_U() const =0;
      
      void generate_momenta(LatticeColourMatrix& P,const RandNum& rand);
    };
    
    namespace MDutils{
      void generate_momenta_su3(LatticeColourMatrix& P,GridParallelRNG& RNG);
    }
    
  }
}
#endif//INTEGRATOR_INCLUDED
