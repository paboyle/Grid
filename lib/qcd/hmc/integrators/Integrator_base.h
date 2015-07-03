//--------------------------------------------------------------------
/*! @file Integrator_base.h
 * @brief Declaration of classes for the abstract Molecular Dynamics integrator
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

class Observer;


/*! @brief Abstract base class for Molecular Dynamics management */

namespace Grid{
  namespace QCD{

    typedef Action<LatticeLorentzColourMatrix>*  ActPtr; // now force the same size as the rest of the code
    typedef std::vector<ActPtr> ActionLevel;
    typedef std::vector<ActionLevel> ActionSet;
    typedef std::vector<Observer*> ObserverList;
    
    
    class Integrator2MN{
      const double lambda = 0.1931833275037836;
      void step (LatticeColourMatrix&, LatticeColourMatrix&, 
		 int, std::vector<int>&);
      
    };
    
    class IntegratorLeapFrog{
      void step (LatticeColourMatrix&, LatticeColourMatrix&,
		 int, std::vector<int>&);
    };
    


    template< class IntegratorPolicy >
    class Integrator{
    private:
      int Nexp;
      int MDsteps;  // number of outer steps
      RealD trajL;  // trajectory length 
      RealD stepsize;

      const std::vector<int> Nrel; // relative steps per level
      const ActionSet as;
      ObserverList observers;
      //      LatticeColourMatrix* const U;  // is shared among all actions - or use a singleton...
      LatticeColourMatrix P;
      
      IntegratorPolicy TheIntegrator;// contains parameters too
      void update_P(int lv,double ep);
      void update_U(double ep);
      
      void register_observers();
      void notify_observers();
      void integrator_step(int level ,std::vector<Integer>& clock);

      
    public:
      Integrator(int Nexp_, int MDsteps_, RealD trajL_,
		 ActionSet& Aset, ObserverList obs):as(Aset), observers(obs){};
      ~Integrator(){}
      void init(LatticeLorentzColourMatrix&,
		GridParallelRNG&);
      double S();
      void integrate(int level);
      LatticeColourMatrix get_U();
    };
    
    namespace MDutils{
      void generate_momenta(LatticeLorentzColourMatrix&,GridParallelRNG&);
      void generate_momenta_su3(LatticeLorentzColourMatrix&,GridParallelRNG&);
    }
    
  }
}
#endif//INTEGRATOR_INCLUDED
