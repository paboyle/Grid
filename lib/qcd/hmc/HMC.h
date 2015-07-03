//--------------------------------------------------------------------
/*! @file HMC.h
 * @brief Declaration of classes for HybridMonteCarlo update
 *
 * @author Guido Cossu
 */
//--------------------------------------------------------------------
#ifndef HMC_INCLUDED
#define HMC_INCLUDED

#include <string>
#include <memory>

namespace Grid{
  namespace QCD{
    
    struct HMCparameters{
      Integer Nsweeps; /* @brief Number of sweeps in this run */
      Integer TotalSweeps; /* @brief If provided, the total number of sweeps */
      Integer ThermalizationSteps;
      Integer StartingConfig;
      Integer SaveInterval; //Setting to 0 does not save configurations
      std::string Filename_prefix; // To save configurations
      
      HMCparameters();
    };
    
 
    class HybridMonteCarlo{
      const HMCparameters Params;
      GridSerialRNG& RNG; 
      // FIXME need the integrator
      
      bool metropolis_test(const RealD DeltaH)const;
      RealD evolve_step(LatticeColourMatrix&)const;
      
      
    public:
      HybridMonteCarlo(GridSerialRNG&);
      ~HybridMonteCarlo(){};

      void evolve(LatticeColourMatrix& Uin)const;

    };
    
  }// QCD
}// Grid


#endif 
