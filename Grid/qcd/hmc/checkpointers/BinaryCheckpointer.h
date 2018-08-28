/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/BinaryCheckpointer.h

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
#ifndef BINARY_CHECKPOINTER
#define BINARY_CHECKPOINTER

#include <iostream>
#include <sstream>
#include <string>

namespace Grid {
namespace QCD {

// Simple checkpointer, only binary file
template <class Impl>
class BinaryHmcCheckpointer : public BaseHmcCheckpointer<Impl> {
 private:
  CheckpointerParameters Params;

 public:
  INHERIT_FIELD_TYPES(Impl);  // Gets the Field type, a Lattice object

  // Extract types from the Field
  typedef typename Field::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;
  typedef typename sobj::DoublePrecision sobj_double;

  BinaryHmcCheckpointer(const CheckpointerParameters &Params_) {
    initialize(Params_);
  }

  void initialize(const CheckpointerParameters &Params_) { Params = Params_; }

  void truncate(std::string file) {
    std::ofstream fout(file, std::ios::out);
    fout.close();
  }

  void TrajectoryComplete(int traj, Field &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG) {

    if ((traj % Params.saveInterval) == 0) {
      std::string config, rng;
      this->build_filenames(traj, Params, config, rng);

      uint32_t nersc_csum;
      uint32_t scidac_csuma;
      uint32_t scidac_csumb;
      
      BinarySimpleUnmunger<sobj_double, sobj> munge;
      truncate(rng);
      BinaryIO::writeRNG(sRNG, pRNG, rng, 0,nersc_csum,scidac_csuma,scidac_csumb);
      truncate(config);

      BinaryIO::writeLatticeObject<vobj, sobj_double>(U, config, munge, 0, Params.format,
						      nersc_csum,scidac_csuma,scidac_csumb);

      std::cout << GridLogMessage << "Written Binary Configuration " << config
                << " checksum " << std::hex 
		<< nersc_csum   <<"/"
		<< scidac_csuma   <<"/"
		<< scidac_csumb 
		<< std::dec << std::endl;
    }

  };

  void CheckpointRestore(int traj, Field &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG) {
    std::string config, rng;
    this->build_filenames(traj, Params, config, rng);
    this->check_filename(rng);
    this->check_filename(config);


    BinarySimpleMunger<sobj_double, sobj> munge;

    uint32_t nersc_csum;
    uint32_t scidac_csuma;
    uint32_t scidac_csumb;
    BinaryIO::readRNG(sRNG, pRNG, rng, 0,nersc_csum,scidac_csuma,scidac_csumb);
    BinaryIO::readLatticeObject<vobj, sobj_double>(U, config, munge, 0, Params.format,
						   nersc_csum,scidac_csuma,scidac_csumb);
    
    std::cout << GridLogMessage << "Read Binary Configuration " << config
              << " checksums " << std::hex << nersc_csum<<"/"<<scidac_csuma<<"/"<<scidac_csumb 
	      << std::dec << std::endl;
  };
};
}
}
#endif
