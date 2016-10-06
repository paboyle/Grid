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

template <class fobj, class sobj>
struct BinarySimpleUnmunger {
  typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

  void operator()(sobj &in, fobj &out, uint32_t &csum) {
    // take word by word and transform accoding to the status
  	fobj_stype* out_buffer = (fobj_stype*)&out;
  	sobj_stype* in_buffer = (sobj_stype*)&in;
  	size_t fobj_words = sizeof(out)/sizeof(fobj_stype);
  	size_t sobj_words = sizeof(in)/sizeof(sobj_stype);
  	assert(fobj_words == sobj_words);

  	for (unsigned int word = 0; word < sobj_words; word++)
			out_buffer[word] = in_buffer[word];  // type conversion on the fly		

			BinaryIO::Uint32Checksum((uint32_t*)&out,sizeof(out),csum);
    
  };

template <class fobj, class sobj>
struct BinarySimpleMunger {
  typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

  void operator()(sobj &out, fobj &in, uint32_t &csum) {
    // take word by word and transform accoding to the status
  	fobj_stype* in_buffer = (fobj_stype*)&in;
  	sobj_stype* out_buffer = (sobj_stype*)&out;
  	size_t fobj_words = sizeof(in)/sizeof(fobj_stype);
  	size_t sobj_words = sizeof(out)/sizeof(sobj_stype);
  	assert(fobj_words == sobj_words);

  	for (unsigned int word = 0; word < sobj_words; word++)
			out_buffer[word] = in_buffer[word];  // type conversion on the fly		

			BinaryIO::Uint32Checksum((uint32_t*)&in,sizeof(in),csum);
    
  };


  // Only for the main field in the hmc
  template <class Impl>
  class BinaryHmcCheckpointer : public HmcObservable<typename Impl::Field> {
   private:
    std::string configStem;
    std::string rngStem;
    int SaveInterval;

   public:
    INHERIT_FIELD_TYPES(Impl);  // The Field is a Lattice object

    typedef typename Field::vector_object vobj;
    typedef typename vobj::scalar_object sobj;
		typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;
		typedef typename sobj::DoublePrecision sobj_double;

    BinaryHmcCheckpointer(std::string cf, std::string rn, int savemodulo, const std::string &format)
        : configStem(cf),
          rngStem(rn),
          SaveInterval(savemodulo){};

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

        // Save always in double precision
        BinarySimpleUnmunger<sobj_double, sobj> munge;
        BinaryIO::writeRNGSerial(sRNG, pRNG, rng, 0);
        BinaryIO::writeObjectParallel<vobj, sobj_double>(U, config, munge, 0, format);
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

      BinarySimpleMunger<sobj_double, sobj> munge;
      BinaryIO::readRNGSerial(sRNG, pRNG, rng, header);
      BinaryIO::readObjectParallel<vobj, sobj_double>(U, config, munge, 0, format);
    };
  };
}
}
#endif
