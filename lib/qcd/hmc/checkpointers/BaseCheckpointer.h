/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/BaseCheckpointer.h

Copyright (C) 2015

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
#ifndef BASE_CHECKPOINTER
#define BASE_CHECKPOINTER

namespace Grid {
namespace QCD {

class CheckpointerParameters : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CheckpointerParameters, 
  	std::string, config_prefix, 
  	std::string, rng_prefix, 
  	int, saveInterval, 
  	std::string, format, );

  CheckpointerParameters(std::string cf = "cfg", std::string rn = "rng",
   		      int savemodulo = 1, const std::string &f = "IEEE64BIG")
      : config_prefix(cf),
        rng_prefix(rn),
        saveInterval(savemodulo),
        format(f){};


  template <class ReaderClass >
  CheckpointerParameters(Reader<ReaderClass> &Reader) {
    read(Reader, "Checkpointer", *this);
  }
 

};

//////////////////////////////////////////////////////////////////////////////
// Base class for checkpointers
template <class Impl>
class BaseHmcCheckpointer : public HmcObservable<typename Impl::Field> {
 public:
  void build_filenames(int traj, CheckpointerParameters &Params,
                       std::string &conf_file, std::string &rng_file) {
    {
      std::ostringstream os;
      os << Params.rng_prefix << "." << traj;
      rng_file = os.str();
    }

    {
      std::ostringstream os;
      os << Params.config_prefix << "." << traj;
      conf_file = os.str();
    }
 	} 

  void check_filename(const std::string &filename){
    std::ifstream f(filename.c_str());
    if(!f.good()){
      std::cout << GridLogError << "Filename " << filename << " not found. Aborting. " << std::endl;
      abort();
    };
  }

  virtual void initialize(const CheckpointerParameters &Params) = 0;

  virtual void CheckpointRestore(int traj, typename Impl::Field &U,
                                 GridSerialRNG &sRNG,
                                 GridParallelRNG &pRNG) = 0;

};  // class BaseHmcCheckpointer
///////////////////////////////////////////////////////////////////////////////
}
}
#endif
