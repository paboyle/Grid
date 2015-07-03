#include <Grid.h>

namespace Grid{
  namespace QCD{

    HMCparameters::HMCparameters(){
	// FIXME fill this constructor  now just default values
	  
	////////////////////////////// Default values
	Nsweeps             = 10;
	TotalSweeps         = 10;
	ThermalizationSteps = 0;
	StartingConfig      = 0;
	SaveInterval        = 1;
	Filename_prefix     = "Conf_";
	/////////////////////////////////
	  
      }


    /////////////////////////////////////////////////////////////////

    HybridMonteCarlo::HybridMonteCarlo(GridSerialRNG& R):RNG(R){
	//FIXME
      }


    void HybridMonteCarlo::evolve(LatticeColourMatrix& Uin)const{
      Real DeltaH;
      Real timer;
      
      // Thermalizations
      for(int iter=1; iter <= Params.ThermalizationSteps; ++iter){
	std::cout << "-- # Thermalization step = "<< iter <<  "\n";
	
	DeltaH = evolve_step(Uin);
	
	std::cout<< "[Timing] Trajectory time (s) : "<< timer/1000.0 << "\n";
	std::cout<< "dH = "<< DeltaH << "\n";
	
	// Update matrix
	//Uin = md_->get_U();  //accept every time
      }

      // Actual updates
      for(int iter=Params.StartingConfig; 
	  iter < Params.Nsweeps+Params.StartingConfig; ++iter){
	std::cout << "-- # Sweep = "<< iter <<  "\n";
		
	DeltaH = evolve_step(Uin);
		
	if(metropolis_test(DeltaH))  {};//Uin = md_->get_U();
	// need sync?	
      }

      
    }


    RealD HybridMonteCarlo::evolve_step(LatticeColourMatrix& Uin)const{
      /*
      md_->init(Uin,RNG);     // set U and initialize P and phi's 
      RealD H0 = md_->S();     // current state            
      std::cout<<"Total H_before = "<< H0 << "\n";
      
      md_->integrator();
      
      RealD H1 = md_->calc_H();     // updated state            
      std::cout<<"Total H_after = "<< H1 << "\n";
      
      return (H1-H0);
      */
      return 0;
    }




    bool HybridMonteCarlo::metropolis_test(const RealD DeltaH)const{
      RealD rn_test;
      RealD prob = std::exp(-DeltaH);
      random(RNG,rn_test);
      
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

  }
}
