/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file:

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

int main(int argc, char **argv)
{
  using namespace Grid;

  Grid_init(&argc, &argv);

  Coordinate latt4  = GridDefaultLatt();
  Coordinate mpi    = GridDefaultMpi();
  Coordinate simd   = GridDefaultSimd(Nd,vComplexD::Nsimd());

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4,simd,mpi);

  GridSerialRNG   sRNG;         sRNG.SeedUniqueString(std::string("The Serial RNG"));
  GridParallelRNG pRNG(UGrid);  pRNG.SeedUniqueString(std::string("The 4D RNG"));

  std::string rngfile("ckpoint_rng.0");
  NerscIO::writeRNGState(sRNG, pRNG, rngfile);
  
  Grid_finalize();
}



