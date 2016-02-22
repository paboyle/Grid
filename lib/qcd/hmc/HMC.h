    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/hmc/HMC.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
//--------------------------------------------------------------------
/*! @file HMC.h
 * @brief Classes for Hybrid Monte Carlo update
 *
 * @author Guido Cossu
 * Time-stamp: <2015-07-30 16:58:26 neo>
 */
//--------------------------------------------------------------------
#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include <string>


namespace Grid{
  namespace QCD{
    

    struct HMCparameters{

      Integer StartTrajectory;
      Integer Trajectories; /* @brief Number of sweeps in this run */
      bool    MetropolisTest;
      Integer NoMetropolisUntil;

      HMCparameters(){
	////////////////////////////// Default values
	MetropolisTest      = true;
	NoMetropolisUntil   = 10;
	StartTrajectory     = 0;
	Trajectories        = 200;
	/////////////////////////////////
      }
    };

    template<class GaugeField> 
    class HmcObservable {
    public:
      virtual void TrajectoryComplete (int traj, GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG & pRNG )=0;
    };

    template<class Gimpl> 
    class PlaquetteLogger : public HmcObservable<typename Gimpl::GaugeField> {
    private:
      std::string Stem;
    public:
      INHERIT_GIMPL_TYPES(Gimpl);
      PlaquetteLogger(std::string cf) {
        Stem  = cf;
      };

      void TrajectoryComplete(int traj, GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG & pRNG )
      {
	  std::string file;   { std::ostringstream os; os << Stem     <<"."<< traj; file = os.str(); }
	  std::ofstream of(file);

	  RealD peri_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(U);
	  RealD peri_rect = WilsonLoops<PeriodicGimplR>::avgRectangle(U);

	  RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(U);
	  RealD impl_rect = WilsonLoops<Gimpl>::avgRectangle(U);

	  of << traj<<" "<< impl_plaq << " " << impl_rect << "  "<< peri_plaq<<" "<<peri_rect<<std::endl;
	  std::cout<< GridLogMessage<< "traj"<<" "<< "plaq " << " " << " rect  " << "  "<< "peri_plaq" <<" "<<"peri_rect"<<std::endl;
	  std::cout<< GridLogMessage<< traj<<" "<< impl_plaq << " " << impl_rect << "  "<< peri_plaq<<" "<<peri_rect<<std::endl;
      }
    };

    //    template <class GaugeField, class Integrator, class Smearer, class Boundary> 
    template <class GaugeField, class IntegratorType>
    class HybridMonteCarlo {
    private:

      const HMCparameters Params;
      
      GridSerialRNG   &sRNG; // Fixme: need a RNG management strategy.
      GridParallelRNG &pRNG; // Fixme: need a RNG management strategy.
      GaugeField      & Ucur;

      IntegratorType &TheIntegrator;
      std::vector<HmcObservable<GaugeField> *> Observables;

      /////////////////////////////////////////////////////////
      // Metropolis step
      /////////////////////////////////////////////////////////
      bool metropolis_test(const RealD DeltaH){

	RealD rn_test;

	RealD prob = std::exp(-DeltaH);

	random(sRNG,rn_test);
      
	std::cout<<GridLogMessage<< "--------------------------------------------\n";
	std::cout<<GridLogMessage<< "dH = "<<DeltaH << "  Random = "<< rn_test <<"\n";
	std::cout<<GridLogMessage<< "Acc. Probability = " << ((prob<1.0)? prob: 1.0)<< "   ";
      
	if((prob >1.0) || (rn_test <= prob)){       // accepted
	  std::cout<<GridLogMessage <<"-- ACCEPTED\n";
	  return true;
	} else {                               // rejected
	  std::cout<<GridLogMessage <<"-- REJECTED\n";
	  return false;
	}

      }

      /////////////////////////////////////////////////////////
      // Evolution
      /////////////////////////////////////////////////////////
      RealD evolve_step(GaugeField& U){

	TheIntegrator.refresh(U,pRNG); // set U and initialize P and phi's 

	RealD H0 = TheIntegrator.S(U); // initial state action  

	std::cout<<GridLogMessage<<"Total H before = "<< H0 << "\n";

	TheIntegrator.integrate(U);
      
	RealD H1 = TheIntegrator.S(U); // updated state action            

	std::cout<<GridLogMessage<<"Total H after = "<< H1 << "\n";

	return (H1-H0);
      }
      
    public:

      /////////////////////////////////////////
      // Constructor
      /////////////////////////////////////////
      HybridMonteCarlo(HMCparameters Pms,  IntegratorType &_Int, GridSerialRNG &_sRNG, GridParallelRNG &_pRNG, GaugeField &_U ) :
        Params(Pms), 
	TheIntegrator(_Int), 
	sRNG(_sRNG),
	pRNG(_pRNG),
	Ucur(_U)
      {
      }
      ~HybridMonteCarlo(){};

      void AddObservable(HmcObservable<GaugeField> *obs) {
	Observables.push_back(obs);
      }

      void evolve(void){

	Real DeltaH;

	GaugeField Ucopy(Ucur._grid);
	
	// Actual updates (evolve a copy Ucopy then copy back eventually)
	for(int traj=Params.StartTrajectory; traj < Params.Trajectories+Params.StartTrajectory; ++traj){

	  std::cout<<GridLogMessage << "-- # Trajectory = "<< traj <<  "\n";
	  Ucopy = Ucur;

	  DeltaH = evolve_step(Ucopy);

	  bool accept = true;
	  if ( traj > Params.NoMetropolisUntil) { 
	    accept = metropolis_test(DeltaH);
	  }
	  
	  if ( accept ) {
	    Ucur = Ucopy;
	  }

	  for(int obs = 0;obs<Observables.size();obs++){
	    Observables[obs]->TrajectoryComplete (traj+1,Ucur,sRNG,pRNG);
	  }

	}
      }
    };
    
  }// QCD
}// Grid


#endif 
