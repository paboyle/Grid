/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/ScidacCheckpointer.h

Copyright (C) 2018

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
#ifndef SCIDAC_CHECKPOINTER
#define SCIDAC_CHECKPOINTER

#ifdef HAVE_LIME

#include <iostream>
#include <sstream>
#include <string>

namespace Grid {
namespace QCD {

// For generic fields
template <class Implementation, class Metadata>
class ScidacHmcCheckpointer : public BaseHmcCheckpointer<Implementation> {
 private:
  CheckpointerParameters Params;
  Metadata MData;

  typedef typename Implementation::Field Field;

 public:
  //INHERIT_GIMPL_TYPES(Implementation);

  ScidacHmcCheckpointer(const CheckpointerParameters &Params_) { initialize(Params_); }
  ScidacHmcCheckpointer(const CheckpointerParameters &Params_, const Metadata& M_):MData(M_) { initialize(Params_); }

  void initialize(const CheckpointerParameters &Params_) {
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

  void TrajectoryComplete(int traj, Field &U, GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {
    if ((traj % Params.saveInterval) == 0) {
      std::string config, rng;
      this->build_filenames(traj, Params, config, rng);
      GridBase *grid = U._grid;
      uint32_t nersc_csum,scidac_csuma,scidac_csumb;
      BinaryIO::writeRNG(sRNG, pRNG, rng, 0,nersc_csum,scidac_csuma,scidac_csumb);
      ScidacWriter _ScidacWriter(grid->IsBoss());
      _ScidacWriter.open(config);
      _ScidacWriter.writeScidacFieldRecord(U, MData);
      _ScidacWriter.close();

      std::cout << GridLogMessage << "Written Scidac Configuration on " << config << std::endl;
    }
  };

  void CheckpointRestore(int traj, Field &U, GridSerialRNG &sRNG,
                         GridParallelRNG &pRNG) {
    std::string config, rng;
    this->build_filenames(traj, Params, config, rng);
    this->check_filename(rng);
    this->check_filename(config);


    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    BinaryIO::readRNG(sRNG, pRNG, rng, 0,nersc_csum,scidac_csuma,scidac_csumb);

    Metadata md_content;
    ScidacReader _ScidacReader;
    _ScidacReader.open(config);
    _ScidacReader.readScidacFieldRecord(U,md_content);  // format from the header
    _ScidacReader.close();

    std::cout << GridLogMessage << "Read Scidac Configuration from " << config
              << " checksum " << std::hex 
	      << nersc_csum<<"/"
	      << scidac_csuma<<"/"
	      << scidac_csumb
	      << std::dec << std::endl;
  };
};
}
}

#endif  // HAVE_LIME
#endif  // ILDG_CHECKPOINTER
