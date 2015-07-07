//--------------------------------------------------------------------
/*! @file HMC.h
 * @brief Classes for Hybrid Monte Carlo update
 *
 * @author Guido Cossu
 * Time-stamp: <2015-07-07 14:58:13 neo>
 */
//--------------------------------------------------------------------
#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include <string>


namespace Grid{
  namespace QCD{
    
    struct HMCparameters{
      Integer Nsweeps; /* @brief Number of sweeps in this run */
      Integer TotalSweeps; /* @brief If provided, the total number of sweeps */
      Integer ThermalizationSteps;
      Integer StartingConfig;
      Integer SaveInterval; //Setting to 0 does not save configurations
      std::string Filename_prefix; // To save configurations and rng seed
      
      HMCparameters();
    };
    
    template <class Algorithm> 
    class HybridMonteCarlo{
      const HMCparameters Params;
      GridSerialRNG sRNG;
      Integrator<Algorithm>& MD;
      
      bool metropolis_test(const RealD DeltaH){
	RealD rn_test;
	RealD prob = std::exp(-DeltaH);
	random(sRNG,rn_test);
      
	std::cout<< "--------------------------------------------\n";
	std::cout<< "dH = "<<DeltaH << "  Random = "<< rn_test 
		 << "\nAcc. Probability = " << ((prob<1.0)? prob: 1.0)<< "   ";
      
	if((prob >1.0) || (rn_test <= prob)){       // accepted
	  std::cout <<"-- ACCEPTED\n";
	  return true;
	} else {                               // rejected
	  std::cout <<"-- REJECTED\n";
	  return false;
	}
      }

      RealD evolve_step(LatticeGaugeField& U){
	MD.init(U); // set U and initialize P and phi's 
	RealD H0 = MD.S(U); // initial state action  
	std::cout<<"Total H before = "<< H0 << "\n";
      
	MD.integrate(U);
      
	RealD H1 = MD.S(U); // updated state action            
	std::cout<<"Total H after = "<< H1 << "\n";
      
	return (H1-H0);
      }
      
    public:
    HybridMonteCarlo(HMCparameters Pms, 
		     Integrator<Algorithm>& MolDyn):
      Params(Pms),MD(MolDyn){
	//FIXME

	// initialize RNGs also with seed
	sRNG.SeedRandomDevice();
      }
      ~HybridMonteCarlo(){};

      void evolve(LatticeGaugeField& Uin){
	Real DeltaH;
	
	// Thermalizations
	for(int iter=1; iter <= Params.ThermalizationSteps; ++iter){
	  std::cout << "-- # Thermalization step = "<< iter <<  "\n";
	
	  DeltaH = evolve_step(Uin);
	  std::cout<< " dH = "<< DeltaH << "\n";
	}

	// Actual updates (evolve a copy Ucopy then copy back eventually)
	LatticeGaugeField Ucopy(Uin._grid);
	for(int iter=Params.StartingConfig; 
	    iter < Params.Nsweeps+Params.StartingConfig; ++iter){
	  std::cout << "-- # Sweep = "<< iter <<  "\n";
	  
	  Ucopy = Uin;
	  DeltaH = evolve_step(Ucopy);
		
	  if(metropolis_test(DeltaH)) Uin = Ucopy;
	  
	  // here save config and RNG seed
	}
      }
    };
    
  }// QCD
}// Grid


#endif 
