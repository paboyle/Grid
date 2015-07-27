#include <Grid.h>

namespace Grid{
  namespace QCD{

    HMCparameters::HMCparameters(){
	// FIXME fill this constructor  now just default values
	  
	////////////////////////////// Default values
	Nsweeps             = 100;
	TotalSweeps         = 20;
	ThermalizationSteps = 20;
	StartingConfig      = 0;
	SaveInterval        = 1;
	Filename_prefix     = "Conf_";
	/////////////////////////////////
	  
      }

  }
}
