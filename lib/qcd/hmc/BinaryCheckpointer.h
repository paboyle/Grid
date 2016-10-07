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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef BINARY_CHECKPOINTER
#define BINARY_CHECKPOINTER

#include <iostream>
#include <sstream>
#include <string>

namespace Grid {
namespace QCD {

// Simple checkpointer, only binary file
template <class Impl>
class BinaryHmcCheckpointer : public HmcObservable<typename Impl::Field> {
 private:
  std::string configStem;
  std::string rngStem;
  int SaveInterval;
  std::string format;

 public:
  INHERIT_FIELD_TYPES(Impl);  // Gets the Field type, a Lattice object

  // Extract types from the Field 
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;
  typedef typename sobj::DoublePrecision sobj_double;

  BinaryHmcCheckpointer(std::string cf, std::string rn, int savemodulo,
                        const std::string &f)
      : configStem(cf), rngStem(rn), SaveInterval(savemodulo), format(f){};

  void truncate(std::string file) {
    std::ofstream fout(file, std::ios::out);
    fout.close();
  }

  void TrajectoryComplete(int traj, Field &U, GridSerialRNG &sRNG,
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

      BinaryIO::BinarySimpleUnmunger<sobj_double, sobj> munge;
      truncate(rng);
      BinaryIO::writeRNGSerial(sRNG, pRNG, rng, 0);
      truncate(config);
      uint32_t csum = BinaryIO::writeObjectParallel<vobj, sobj_double>(
          U, config, munge, 0, format);

      std::cout << GridLogMessage << "Written Binary Configuration " << config
                << " checksum " << std::hex << csum << std::dec << std::endl;
    }
  };

  void CheckpointRestore(int traj, Field &U, GridSerialRNG &sRNG,
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

    BinaryIO::BinarySimpleMunger<sobj_double, sobj> munge;
    BinaryIO::readRNGSerial(sRNG, pRNG, rng, 0);
    uint32_t csum = BinaryIO::readObjectParallel<vobj, sobj_double>(
        U, config, munge, 0, format);

    std::cout << GridLogMessage << "Read Binary Configuration " << config
              << " checksum " << std::hex << csum << std::dec << std::endl;
  };
};
}
}
#endif
