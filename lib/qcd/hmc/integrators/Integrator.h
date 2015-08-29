//--------------------------------------------------------------------
/*! @file Integrator.h
 * @brief Classes for the Molecular Dynamics integrator
 *
 * @author Guido Cossu
 * Time-stamp: <2015-07-30 16:21:29 neo>
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

class Observer;

#include <memory>

namespace Grid{
  namespace QCD{

    typedef Action<LatticeGaugeField>*  ActPtr; // now force the same colours as the rest of the code
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
      void generate_momenta(LatticeGaugeField&,GridParallelRNG&);
      void generate_momenta_su3(LatticeGaugeField&,GridParallelRNG&);
    }

    /*! @brief Class for Molecular Dynamics management */   
    template< class IntegratorAlgorithm >
    class Integrator{
    private:

      int levels;
      std::vector<double> t_P;
      double t_U;

      IntegratorParameters Params;
      const ActionSet as;
      std::unique_ptr<LatticeGaugeField> P;
      GridParallelRNG pRNG;

      //ObserverList observers; // not yet
     
      IntegratorAlgorithm TheIntegrator;

      void register_observers();

      void notify_observers();

      void update_P(LatticeGaugeField&U, int level,double ep){
	t_P[level]+=ep;
	for(int a=0; a<as[level].actions.size(); ++a){
	  LatticeGaugeField force(U._grid);
	  as[level].actions.at(a)->deriv(U,force);
	  *P -= force*ep;
	}
	std::cout<<GridLogMessage;
	for(int l=0; l<level;++l) std::cout<<"   ";	    
	std::cout<<"["<<level<<"] P " << " dt "<< ep <<" : t_P "<< t_P[level] <<std::endl;
      }

      void update_U(LatticeGaugeField&U, double ep){
	//rewrite exponential to deal automatically  with the lorentz index?
	LatticeColourMatrix Umu(U._grid);
	LatticeColourMatrix Pmu(U._grid);
	for (int mu = 0; mu < Nd; mu++){
	  Umu=PeekIndex<LorentzIndex>(U, mu);
	  Pmu=PeekIndex<LorentzIndex>(*P, mu);
	  Umu = expMat(Pmu, ep, Params.Nexp)*Umu;
	  PokeIndex<LorentzIndex>(U, Umu, mu);
	}
	t_U+=ep;
	int fl = levels-1;
	std::cout<<GridLogMessage<<"   ";
	for(int l=0; l<fl;++l) std::cout<<"   ";	    
	std::cout<<"["<<fl<<"] U " << " dt "<< ep <<" : t_U "<< t_U <<std::endl;
      }
      
      friend void IntegratorAlgorithm::step (LatticeGaugeField& U, 
					     int level, std::vector<int>& clock,
					     Integrator<IntegratorAlgorithm>* Integ);
    public:

      Integrator(GridBase* grid, 
		 IntegratorParameters Par,
		 ActionSet& Aset):
          Params(Par),
    	  as(Aset),
	  P(new LatticeGaugeField(grid)),
	  pRNG(grid),
	  levels(Aset.size())
      {
	
	std::vector<int> seeds({1,2,3,4,5});
	pRNG.SeedFixedIntegers(seeds);

	t_P.resize(levels,0.0);
	t_U=0.0;

      };
      
      ~Integrator(){}

      //Initialization of momenta and actions
      void init(LatticeGaugeField& U){
	std::cout<<GridLogMessage<< "Integrator init\n";
	MDutils::generate_momenta(*P,pRNG);
	for(int level=0; level< as.size(); ++level){
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    as[level].actions.at(actionID)->init(U, pRNG);
	  }
	}
      }

      // Calculate action
      RealD S(LatticeGaugeField& U){
	LatticeComplex Hloc(U._grid);
	Hloc = zero;
	// Momenta
	for (int mu=0; mu <Nd; mu++){
	  LatticeColourMatrix Pmu = peekLorentz(*P, mu);
	  Hloc -= trace(Pmu*Pmu);
	}
	Complex Hsum = sum(Hloc);
	
	RealD H = Hsum.real();
	RealD Hterm;
	std::cout<<GridLogMessage << "Momentum action H_p = "<< H << "\n";

	// Actions
	for(int level=0; level<as.size(); ++level){
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    Hterm = as[level].actions.at(actionID)->S(U);
	    std::cout<<GridLogMessage << "Level "<<level<<" term "<<actionID<<" H = "<<Hterm<<std::endl;
	    H += Hterm;
	  }
	}
	std::cout<<GridLogMessage << "Total action H = "<< H << "\n";
	
	return H;
      }

      void integrate(LatticeGaugeField& U){
	std::vector<int> clock;
	clock.resize(as.size(),0);
	for(int step=0; step< Params.MDsteps; ++step){   // MD step
	  TheIntegrator.step(U,0,clock, (this));
	}
	for(int level=0; level<as.size(); ++level){
	  assert(fabs(t_U - t_P[level])<1.0e-3);
	  std::cout<<GridLogMessage<<" times["<<level<<"]= "<<t_P[level]<< " " << t_U <<std::endl;
	}	
      }
    };
    
  }
}
#endif//INTEGRATOR_INCLUDED
