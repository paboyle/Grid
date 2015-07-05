//--------------------------------------------------------------------
/*! @file Integrator_base.h
 * @brief Declaration of classes for the abstract Molecular Dynamics integrator
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

#include <memory>

class Observer;


/*! @brief Abstract base class for Molecular Dynamics management */

namespace Grid{
  namespace QCD{

    typedef Action<LatticeLorentzColourMatrix>*  ActPtr; // now force the same size as the rest of the code
    typedef std::vector<ActPtr> ActionLevel;
    typedef std::vector<ActionLevel> ActionSet;
    typedef std::vector<Observer*> ObserverList;
    
    class LeapFrog;
   

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

    
    


    template< class IntegratorPolicy >
    class Integrator{
    private:
      IntegratorParameters Params;
      const ActionSet as;
      const std::vector<int> Nrel; //relative step size per level
      //ObserverList observers; // not yet
      std::unique_ptr<LatticeLorentzColourMatrix> P;
      
      IntegratorPolicy TheIntegrator;// contains parameters too


      void update_P(LatticeLorentzColourMatrix&U, int level,double ep){
	for(int a=0; a<as[level].size(); ++a){
	  LatticeLorentzColourMatrix force(U._grid);
	  as[level].at(a)->deriv(U,force);
	  *P -= force*ep;
	}
      }


      void update_U(LatticeLorentzColourMatrix&U, double ep){
	//rewrite exponential to deal with the lorentz index?
	LatticeColourMatrix Umu(U._grid);
	LatticeColourMatrix Pmu(U._grid);
	for (int mu = 0; mu < Nd; mu++){
	  Umu=peekLorentz(U, mu);
	  Pmu=peekLorentz(*P, mu);
	  std::cout << "U norm ["<<mu<<"] = "<< norm2(Umu) << "  : Pmu "<< norm2(Pmu)<<"\n";
	  Umu = expMat(Pmu, Complex(ep, 0.0))*Umu;
	  pokeLorentz(U, Umu, mu);
	}

      }
      
      void register_observers();
      void notify_observers();
      
      friend void IntegratorPolicy::step (LatticeLorentzColourMatrix& U, 
				     int level, std::vector<int>& clock,
				     Integrator<LeapFrog>* Integ);
    public:
      Integrator(IntegratorParameters Par,
		 ActionSet& Aset, std::vector<int> Nrel_):
	Params(Par),as(Aset),Nrel(Nrel_){
	assert(as.size() == Nrel.size());
      };
      
      ~Integrator(){}


      //Initialization of momenta and actions
      void init(LatticeLorentzColourMatrix& U,
		GridParallelRNG& pRNG){
	std::cout<< "Integrator init\n";
	if (!P){
	  std::unique_ptr<LatticeLorentzColourMatrix> Pnew(new LatticeLorentzColourMatrix(U._grid));
	  P = std::move(Pnew);

	}
	MDutils::generate_momenta(*P,pRNG);
	
       
	
	for(int level=0; level< as.size(); ++level){
	  for(int actionID=0; actionID<as.at(level).size(); ++actionID){
	    as[level].at(actionID)->init(U, pRNG);
	  }
	}
      }

      

      RealD S(LatticeLorentzColourMatrix& U){
	LatticeComplex Hloc(U._grid);
	Hloc = zero;
	// Momenta
	for (int mu=0; mu <Nd; mu++){
	  LatticeColourMatrix Pmu = peekLorentz(*P, mu);
	  Hloc -= trace(Pmu*adj(Pmu));
	}
	Complex Hsum = sum(Hloc);
	
	RealD H = Hsum.real();

	std::cout << "H_p = "<< H << "\n";

	// Actions
	for(int level=0; level<as.size(); ++level)
	  for(int actionID=0; actionID<as.at(level).size(); ++actionID)
	    H += as[level].at(actionID)->S(U);
	
	return H;
      }
      
      

      void integrate(LatticeLorentzColourMatrix& U, int level){
	std::vector<int> clock;
	clock.resize(as.size(),0);
	for(int step=0; step< Params.MDsteps; ++step)   // MD step
	  TheIntegrator.step(U,0,clock, (this));
      }
    };
    

    class MinimumNorm2{
      const double lambda = 0.1931833275037836;
    public:
      void step (LatticeLorentzColourMatrix& U, int level, std::vector<int>& clock);
      
    };
    
    class LeapFrog{
    public:
      void step (LatticeLorentzColourMatrix& U, 
		 int level, std::vector<int>& clock,
		 Integrator<LeapFrog>* Integ){
	// cl  : current level
	// fl  : final level
	// eps : current step size
	
	int fl = Integ->as.size() -1;
	double eps = Integ->Params.stepsize;
	
	// Get current level step size
	for(int l=0; l<=level; ++l) eps/= Integ->Nrel[l];
	
	int fin = 1;
	for(int l=0; l<=level; ++l) fin*= Integ->Nrel[l];
	fin = 2*Integ->Params.MDsteps*fin - 1;
	
	for(int e=0; e<Integ->Nrel[level]; ++e){
	  
	  if(clock[level] == 0){    // initial half step
	    Integ->update_P(U, level,eps/2);
	    ++clock[level];
	    for(int l=0; l<level;++l) std::cout<<"   ";
	    std::cout<<"P "<< 0.5*clock[level] <<std::endl;
	  }
	  if(level == fl){          // lowest level
	    Integ->update_U(U, eps);
	    for(int l=0; l<level;++l) std::cout<<"   ";
	    std::cout<<"U "<< 0.5*(clock[level]+1) <<std::endl;
	  }else{                 // recursive function call
	    step(U, level+1,clock, Integ);
	  }
	  if(clock[level] == fin){  // final half step
	    Integ->update_P(U, level,eps/2);
	    
	    ++clock[level];
	    for(int l=0; l<level;++l) std::cout<<"   ";
	    std::cout<<"P "<< 0.5*clock[level] <<std::endl;
	  }else{                  // bulk step
	    Integ->update_P(U, level,eps);
	    
	    clock[level]+=2;
	    for(int l=0; l<level;++l) std::cout<<"   ";
	    std::cout<<"P "<< 0.5*clock[level] <<std::endl;
	  }
	}




      }
    };




    
  }
}
#endif//INTEGRATOR_INCLUDED
