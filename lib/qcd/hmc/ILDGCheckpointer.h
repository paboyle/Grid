/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/ILDGCheckpointer.h

Copyright (C) 2015

Author: Guido Cossu

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef ILDG_CHECKPOINTER
#define ILDG_CHECKPOINTER

#include <iostream>
#include <sstream>
#include <string>

namespace Grid {
namespace QCD {

// Only for Gauge fields
template <class Implementation>
class ILDGHmcCheckpointer
    : public HmcObservable<typename Implementation::GaugeField> {
 private:
  std::string configStem;
  std::string rngStem;
  int SaveInterval;
  std::string format;

 public:
  INHERIT_GIMPL_TYPES(Implementation);  //

  ILDGHmcCheckpointer(std::string cf, std::string rn, int savemodulo,
                       std::string form = "IEEE64BIG") {
    configStem = cf;
    rngStem = rn;
    SaveInterval = savemodulo;
    format = form;

    // check here that the format is valid
    int ieee32big = (format == std::string("IEEE32BIG"));
    int ieee32    = (format == std::string("IEEE32"));
    int ieee64big = (format == std::string("IEEE64BIG"));
    int ieee64    = (format == std::string("IEEE64"));

    if (!(ieee64big ^ ieee32 ^ ieee32big ^ ieee64)) {
      std::cout << GridLogMessage << "Invalid format: " << format << std::endl;
      exit(0);
    }
  };

  void TrajectoryComplete(int traj, GaugeField &U, GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    if ((traj % SaveInterval) == 0) {
      std::string rng;
      {
        std::ostringstream os;
        os << rngStem << "." << traj;
        rng = os.str();
      }
      std::string config;
      {
        std::ostringstream os;
        os << configStem << "." << traj;
        config = os.str();
      }

      ILDGIO IO(config, ILDGwrite);
      BinaryIO::writeRNGSerial(sRNG, pRNG, rng, 0);
      uint32_t csum  = IO.writeConfiguration(U, format);

      std::cout << GridLogMessage << "Written ILDG Configuration on " << config
                << " checksum " << std::hex << csum << std::dec << std::endl;
    }
  };

  void CheckpointRestore(int traj, GaugeField &U, GridSerialRNG &sRNG,
                         GridParallelRNG &pRNG) {
    std::string rng;
    {
      std::ostringstream os;
      os << rngStem << "." << traj;
      rng = os.str();
    }
    std::string config;
    {
      std::ostringstream os;
      os << configStem << "." << traj;
      config = os.str();
    }

    ILDGIO IO(config, ILDGread);
    BinaryIO::readRNGSerial(sRNG, pRNG, rng, 0);
    uint32_t csum = IO.readConfiguration(U);// format from the header

    std::cout << GridLogMessage << "Read ILDG Configuration from " << config
              << " checksum " << std::hex << csum << std::dec << std::endl;
  };
};
}
}
#endif
