//--------------------------------------------------------------------
/*! @file Integrator.h
 * @brief Declaration of classes for the Molecular Dynamics integrator
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

class Observer;

#include <memory>

namespace Grid{
  namespace QCD{

    typedef Action<LatticeLorentzColourMatrix>*  ActPtr; // now force the same colours as the rest of the code
    struct ActionLevel{
      int multiplier;
    public:
      std::vector<ActPtr> actions;
      explicit ActionLevel(int mul = 1):multiplier(mul){assert (mul > 0);};
      void push_back(ActPtr ptr){
	actions.push_back(ptr);
      }
    };
    typedef std::vector<ActionLevel> ActionSet;
    typedef std::vector<Observer*> ObserverList;
    
    struct IntegratorParameters{
      int Nexp;
      int MDsteps;  // number of outer steps
      RealD trajL;  // trajectory length 
      RealD stepsize;

      IntegratorParameters(int Nexp_,
			   int MDsteps_, 
			   RealD trajL_):
      Nexp(Nexp_),MDsteps(MDsteps_),trajL(trajL_),stepsize(trajL/MDsteps){};
    };


    namespace MDutils{
      void generate_momenta(LatticeLorentzColourMatrix&,GridParallelRNG&);
      void generate_momenta_su3(LatticeLorentzColourMatrix&,GridParallelRNG&);
    }

    /*! @brief Class for Molecular Dynamics management */   
    template< class IntegratorAlgorithm >
    class Integrator{
    private:
      IntegratorParameters Params;
      const ActionSet as;
      std::unique_ptr<LatticeLorentzColourMatrix> P;
      GridParallelRNG pRNG;
      //ObserverList observers; // not yet
     
      IntegratorAlgorithm TheIntegrator;

      void register_observers();
      void notify_observers();

      void update_P(LatticeLorentzColourMatrix&U, int level,double ep){
	for(int a=0; a<as[level].actions.size(); ++a){
	  LatticeLorentzColourMatrix force(U._grid);
	  as[level].actions.at(a)->deriv(U,force);
	  *P -= force*ep;
	}
      }


      void update_U(LatticeLorentzColourMatrix&U, double ep){
	//rewrite exponential to deal automatically  with the lorentz index?
	LatticeColourMatrix Umu(U._grid);
	LatticeColourMatrix Pmu(U._grid);
	for (int mu = 0; mu < Nd; mu++){
	  Umu=peekLorentz(U, mu);
	  Pmu=peekLorentz(*P, mu);
	  Umu = expMat(Pmu, ep, Params.Nexp)*Umu;
	  pokeLorentz(U, Umu, mu);
	}

      }
      

      
      friend void IntegratorAlgorithm::step (LatticeLorentzColourMatrix& U, 
					     int level, std::vector<int>& clock,
					     Integrator<IntegratorAlgorithm>* Integ);
    public:
    Integrator(GridBase* grid, IntegratorParameters Par,
		 ActionSet& Aset):
      Params(Par),as(Aset),P(new LatticeLorentzColourMatrix(grid)),pRNG(grid){
	pRNG.SeedRandomDevice();
      };
      
      ~Integrator(){}


      //Initialization of momenta and actions
      void init(LatticeLorentzColourMatrix& U){
	std::cout<< "Integrator init\n";

	MDutils::generate_momenta(*P,pRNG);
	for(int level=0; level< as.size(); ++level){
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    as[level].actions.at(actionID)->init(U, pRNG);
	  }
	}
      }

      
      // Calculate action
      RealD S(LatticeLorentzColourMatrix& U){
	LatticeComplex Hloc(U._grid);
	Hloc = zero;
	// Momenta
	for (int mu=0; mu <Nd; mu++){
	  LatticeColourMatrix Pmu = peekLorentz(*P, mu);
	  Hloc -= trace(Pmu*Pmu);
	}
	Complex Hsum = sum(Hloc);
	
	RealD H = Hsum.real();

	std::cout << "H_p = "<< H << "\n";

	// Actions
	for(int level=0; level<as.size(); ++level)
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID)
	    H += as[level].actions.at(actionID)->S(U);
	
	return H;
      }

      void integrate(LatticeLorentzColourMatrix& U){
	std::vector<int> clock;
	clock.resize(as.size(),0);
	for(int step=0; step< Params.MDsteps; ++step)   // MD step
	  TheIntegrator.step(U,0,clock, (this));
      }
    };
    





    
  }
}
#endif//INTEGRATOR_INCLUDED
