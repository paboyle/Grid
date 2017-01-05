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
  GRID_SERIALIZABLE_CLASS_MEMBERS(CheckpointerParameters, std::string,
                                  configStem, std::string, rngStem, int,
                                  SaveInterval, std::string, format, );

  CheckpointerParameters(std::string cf = "cfg", std::string rn = "rng",
                         int savemodulo = 1, const std::string &f = "IEEE64BIG")
      : configStem(cf), rngStem(rn), SaveInterval(savemodulo), format(f){};
};

//////////////////////////////////////////////////////////////////////////////
// Base class for checkpointers
template <class Impl>
class BaseHmcCheckpointer : public HmcObservable<typename Impl::Field> {
 public:
  virtual void initialize(CheckpointerParameters &Params) = 0;

  virtual void CheckpointRestore(int traj, typename Impl::Field &U,
                                 GridSerialRNG &sRNG,
                                 GridParallelRNG &pRNG) = 0;

};  // class BaseHmcCheckpointer
///////////////////////////////////////////////////////////////////////////////
}
}
#endif
