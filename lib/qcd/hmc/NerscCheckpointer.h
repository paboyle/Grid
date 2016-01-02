    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/hmc/NerscCheckpointer.h

    Copyright (C) 2015

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
