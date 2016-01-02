#ifndef NERSC_CHECKPOINTER
#define NERSC_CHECKPOINTER

#include <string>
#include <iostream>
#include <sstream>


namespace Grid{
  namespace QCD{
    

    template<class Gimpl> 
    class NerscHmcCheckpointer : public HmcObservable<typename Gimpl::GaugeField> {
    private:
      std::string configStem;
      std::string rngStem;
      int SaveInterval;
    public:
      INHERIT_GIMPL_TYPES(Gimpl);

      NerscHmcCheckpointer(std::string cf, std::string rn,int savemodulo) {
        configStem  = cf;
        rngStem     = rn;
        SaveInterval= savemodulo;
      };

      void TrajectoryComplete(int traj, GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG & pRNG )
      {
	if ( (traj % SaveInterval)== 0 ) {
	  std::string rng;   { std::ostringstream os; os << rngStem     <<"."<< traj; rng = os.str(); }
	  std::string config;{ std::ostringstream os; os << configStem  <<"."<< traj; config = os.str();}

	  int precision32=1;
	  int tworow     =0;
	  NerscIO::writeRNGState(sRNG,pRNG,rng);
	  NerscIO::writeConfiguration(U,config,tworow,precision32);
	}
      };

      void CheckpointRestore(int traj, GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG & pRNG ){

	std::string rng;   { std::ostringstream os; os << rngStem     <<"."<< traj; rng = os.str(); }
	std::string config;{ std::ostringstream os; os << configStem  <<"."<< traj; config = os.str();}

	NerscField header;
	NerscIO::readRNGState(sRNG,pRNG,header,rng);
	NerscIO::readConfiguration(U,header,config);
      };

    };
}}
#endif
