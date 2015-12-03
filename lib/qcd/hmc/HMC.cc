#include <Grid.h>

namespace Grid{
  namespace QCD{

    HMCparameters::HMCparameters(){
	// FIXME fill this constructor  now just default values
	  
	////////////////////////////// Default values
	Nsweeps             = 200;
	TotalSweeps         = 240;
	ThermalizationSteps = 40;
	StartingConfig      = 0;
	SaveInterval        = 1;
	Filename_prefix     = "Conf_";
	/////////////////////////////////
	  
      }

  }
}
