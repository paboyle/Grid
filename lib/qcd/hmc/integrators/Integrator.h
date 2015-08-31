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

//class Observer;

#include <memory>

namespace Grid{
  namespace QCD{

    struct IntegratorParameters{

      int Nexp;
      int MDsteps;  // number of outer steps
      RealD trajL;  // trajectory length 
      RealD stepsize;

      IntegratorParameters(int MDsteps_, 
			   RealD trajL_=1.0,
			   int Nexp_=12):
        Nexp(Nexp_),
	MDsteps(MDsteps_),
	trajL(trajL_),
	stepsize(trajL/MDsteps)
        {
	  // empty body constructor
	};

    };

    /*! @brief Class for Molecular Dynamics management */   
    template<class GaugeField>
    class Integrator {

    protected:

      typedef IntegratorParameters ParameterType;

      IntegratorParameters Params;

      const ActionSet<GaugeField> as;

      int levels;              //
      double t_U;              // Track time passing on each level and for U and for P
      std::vector<double> t_P; //

      GaugeField P;

      // Should match any legal (SU(n)) gauge field
      // Need to use this template to match Ncol to pass to SU<N> class
      template<int Ncol,class vec> void generate_momenta(Lattice< iVector< iScalar< iMatrix<vec,Ncol> >, Nd> > & P,GridParallelRNG& pRNG){
	typedef Lattice< iScalar< iScalar< iMatrix<vec,Ncol> > > > GaugeLinkField;
	GaugeLinkField Pmu(P._grid);
	Pmu = zero;
	for(int mu=0;mu<Nd;mu++){
	  SU<Ncol>::GaussianLieAlgebraMatrix(pRNG, Pmu);
	  PokeIndex<LorentzIndex>(P, Pmu, mu);
	}
      }


      //ObserverList observers; // not yet
      //      typedef std::vector<Observer*> ObserverList;
      //      void register_observers();
      //      void notify_observers();


      void update_P(GaugeField&U, int level,double ep){
	t_P[level]+=ep;
	update_P(P,U,level,ep);

	std::cout<<GridLogMessage;
	for(int l=0; l<level;++l) std::cout<<"   ";	    
	std::cout<<"["<<level<<"] P " << " dt "<< ep <<" : t_P "<< t_P[level] <<std::endl;
      }
      void update_P(GaugeField &Mom,GaugeField&U, int level,double ep){
	for(int a=0; a<as[level].actions.size(); ++a){
	  GaugeField force(U._grid);
	  as[level].actions.at(a)->deriv(U,force);
	  Mom = Mom - force*ep;
	}
      }

      void update_U(GaugeField&U, double ep){
	update_U(P,U,ep);

	t_U+=ep;
	int fl = levels-1;
	std::cout<<GridLogMessage<<"   ";
	for(int l=0; l<fl;++l) std::cout<<"   ";	    
	std::cout<<"["<<fl<<"] U " << " dt "<< ep <<" : t_U "<< t_U <<std::endl;

      }
      void update_U(GaugeField &Mom, GaugeField&U, double ep){
	//rewrite exponential to deal automatically  with the lorentz index?
	//	GaugeLinkField Umu(U._grid);
	//	GaugeLinkField Pmu(U._grid);
	for (int mu = 0; mu < Nd; mu++){
	  auto Umu=PeekIndex<LorentzIndex>(U, mu);
	  auto Pmu=PeekIndex<LorentzIndex>(Mom, mu);
	  Umu = expMat(Pmu, ep, Params.Nexp)*Umu;
	  PokeIndex<LorentzIndex>(U, Umu, mu);
	}
      }
      
      /*
	friend void Algorithm::step (GaugeField& U, 
				   int level, 
				   std::vector<int>& clock,
				   Integrator<GaugeField,Algorithm>* Integ);
      */

      virtual void step (GaugeField& U,int level, std::vector<int>& clock)=0;

    public:

      Integrator(GridBase* grid, 
		 IntegratorParameters Par,
		 ActionSet<GaugeField> & Aset):
          Params(Par),
    	  as(Aset),
	  P(grid),
	  levels(Aset.size())
      {
	t_P.resize(levels,0.0);
	t_U=0.0;
      };
      
      virtual ~Integrator(){}

      //Initialization of momenta and actions
      void refresh(GaugeField& U,GridParallelRNG &pRNG){
	std::cout<<GridLogMessage<< "Integrator refresh\n";
	generate_momenta(P,pRNG);
	for(int level=0; level< as.size(); ++level){
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    as[level].actions.at(actionID)->refresh(U, pRNG);
	  }
	}
      }

      // Calculate action
      RealD S(GaugeField& U){

	LatticeComplex Hloc(U._grid);	Hloc = zero;
	// Momenta
	for (int mu=0; mu <Nd; mu++){
	  auto Pmu = PeekIndex<LorentzIndex>(P, mu);
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

      void integrate(GaugeField& U){

	std::vector<int> clock;

	clock.resize(as.size(),0);

	// All the clock stuff is removed if we pass first, last to the step down the way
	for(int step=0; step< Params.MDsteps; ++step){   // MD step
	  int first_step = (step==0);
	  int  last_step = (step==Params.MDsteps-1);
	  this->step(U,0,clock);
	}

	// Check the clocks all match
	for(int level=0; level<as.size(); ++level){
	  assert(fabs(t_U - t_P[level])<1.0e-6); // must be the same
	  std::cout<<GridLogMessage<<" times["<<level<<"]= "<<t_P[level]<< " " << t_U <<std::endl;
	}	

      }
    };
    
  }
}
#endif//INTEGRATOR_INCLUDED
