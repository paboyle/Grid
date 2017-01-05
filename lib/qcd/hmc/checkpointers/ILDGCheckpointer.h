/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/ILDGCheckpointer.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifdef HAVE_LIME

#include <iostream>
#include <sstream>
#include <string>

namespace Grid {
namespace QCD {

// Only for Gauge fields
template <class Implementation>
class ILDGHmcCheckpointer
    : public BaseHmcCheckpointer<Implementation> {
 private:
 	CheckpointerParameters Params;
/*
  std::string configStem;
  std::string rngStem;
  int SaveInterval;
  std::string format;
*/

 public:
  INHERIT_GIMPL_TYPES(Implementation);

  ILDGHmcCheckpointer(CheckpointerParameters &Params_) { initialize(Params_); }

  void initialize(CheckpointerParameters &Params_) {
    Params = Params_;

    // check here that the format is valid
    int ieee32big = (Params.format == std::string("IEEE32BIG"));
    int ieee32    = (Params.format == std::string("IEEE32"));
    int ieee64big = (Params.format == std::string("IEEE64BIG"));
    int ieee64    = (Params.format == std::string("IEEE64"));

    if (!(ieee64big || ieee32 || ieee32big || ieee64)) {
      std::cout << GridLogError << "Unrecognized file format " << Params.format
                << std::endl;
      std::cout << GridLogError
                << "Allowed: IEEE32BIG | IEEE32 | IEEE64BIG | IEEE64"
                << std::endl;

      exit(1);
    }
  }

  void TrajectoryComplete(int traj, GaugeField &U, GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    if ((traj % Params.SaveInterval) == 0) {
      std::string rng;
      {
        std::ostringstream os;
        os << Params.rngStem << "." << traj;
        rng = os.str();
      }
      std::string config;
      {
        std::ostringstream os;
        os << Params.configStem << "." << traj;
        config = os.str();
      }

      ILDGIO IO(config, ILDGwrite);
      BinaryIO::writeRNGSerial(sRNG, pRNG, rng, 0);
      uint32_t csum  = IO.writeConfiguration(U, Params.format);

      std::cout << GridLogMessage << "Written ILDG Configuration on " << config
                << " checksum " << std::hex << csum << std::dec << std::endl;
    }
  };

  void CheckpointRestore(int traj, GaugeField &U, GridSerialRNG &sRNG,
                         GridParallelRNG &pRNG) {
    std::string rng;
    {
      std::ostringstream os;
      os << Params.rngStem << "." << traj;
      rng = os.str();
    }
    std::string config;
    {
      std::ostringstream os;
      os << Params.configStem << "." << traj;
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

#endif // HAVE_LIME
#endif // ILDG_CHECKPOINTER
