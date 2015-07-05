#include <Grid.h>

namespace Grid{
  namespace QCD{

    HMCparameters::HMCparameters(){
	// FIXME fill this constructor  now just default values
	  
	////////////////////////////// Default values
	Nsweeps             = 10;
	TotalSweeps         = 10;
	ThermalizationSteps = 1;
	StartingConfig      = 0;
	SaveInterval        = 1;
	Filename_prefix     = "Conf_";
	/////////////////////////////////
	  
      }


   

  }
}
