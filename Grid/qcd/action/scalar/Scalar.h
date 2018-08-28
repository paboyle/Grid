/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/Scalar.h

Copyright (C) 2017

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
#ifndef GRID_QCD_SCALAR_H
#define GRID_QCD_SCALAR_H

#include <Grid/qcd/action/scalar/ScalarImpl.h>
#include <Grid/qcd/action/scalar/ScalarAction.h>
#include <Grid/qcd/action/scalar/ScalarInteractionAction.h>

namespace Grid {
namespace QCD {

  typedef ScalarAction<ScalarImplR>                 ScalarActionR;
  typedef ScalarAction<ScalarImplF>                 ScalarActionF;
  typedef ScalarAction<ScalarImplD>                 ScalarActionD;

  template <int Colours, int Dimensions> using ScalarAdjActionR = ScalarInteractionAction<ScalarNxNAdjImplR<Colours>, Dimensions>;
  template <int Colours, int Dimensions> using ScalarAdjActionF = ScalarInteractionAction<ScalarNxNAdjImplF<Colours>, Dimensions>;
  template <int Colours, int Dimensions> using ScalarAdjActionD = ScalarInteractionAction<ScalarNxNAdjImplD<Colours>, Dimensions>;
  
}
}

#endif  // GRID_QCD_SCALAR_H
