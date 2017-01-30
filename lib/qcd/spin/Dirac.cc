/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: lib/qcd/spin/Dirac.cc

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#include <Grid.h>

namespace Grid {
namespace QCD {

#include "GammaMulTable.h"

const std::array<const char *, Gamma::nGamma> Gamma::name = {{
  "-Gamma5      ",
  "Gamma5       ",
  "-GammaT      ",
  "GammaT       ",
  "-GammaTGamma5",
  "GammaTGamma5 ",
  "-GammaX      ",
  "GammaX       ",
  "-GammaXGamma5",
  "GammaXGamma5 ",
  "-GammaY      ",
  "GammaY       ",
  "-GammaYGamma5",
  "GammaYGamma5 ",
  "-GammaZ      ",
  "GammaZ       ",
  "-GammaZGamma5",
  "GammaZGamma5 ",
  "-Identity    ",
  "Identity     ",
  "-SigmaXT     ",
  "SigmaXT      ",
  "-SigmaXY     ",
  "SigmaXY      ",
  "-SigmaXZ     ",
  "SigmaXZ      ",
  "-SigmaYT     ",
  "SigmaYT      ",
  "-SigmaYZ     ",
  "SigmaYZ      ",
  "-SigmaZT     ",
  "SigmaZT      "}};

}}
