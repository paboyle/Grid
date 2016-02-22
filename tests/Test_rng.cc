    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_rng.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
     
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});

  GridSerialRNG             sRNG;   sRNG.SeedRandomDevice();
  GridSerialRNG            fsRNG;  fsRNG.SeedFixedIntegers(seeds);

  GridParallelRNG           pRNG(&Grid);   pRNG.SeedRandomDevice();
  GridParallelRNG          fpRNG(&Grid);  fpRNG.SeedFixedIntegers(seeds);

  SpinMatrix rnd  ; 
  random(sRNG,rnd);
  std::cout<<GridLogMessage<<"Random Spin Matrix (random_device)\n"<< rnd<<std::endl;

  random(fsRNG,rnd);
  std::cout<<GridLogMessage<<"Random Spin Matrix (fixed seed)\n"<< rnd<<std::endl;

  SpinVector rv; 
  random(sRNG,rv);
  std::cout<<GridLogMessage<<"Random Spin Vector (random device)\n"<< rv<<std::endl;

  random(fsRNG,rv);
  std::cout<<GridLogMessage<<"Random Spin Vector (fixed seed)\n"<< rv<<std::endl;

  gaussian(fsRNG,rv);
  std::cout<<GridLogMessage<<"Gaussian Spin Vector (fixed seed)\n"<< rv<<std::endl;

  LatticeColourVector lcv(&Grid);
  random(pRNG,lcv);
  std::cout<<GridLogMessage<<"Random Lattice Colour Vector (random device)\n"<< lcv<<std::endl;

  random(fpRNG,lcv);
  std::cout<<GridLogMessage<<"Random Lattice Colour Vector (fixed seed)\n"<< lcv<<std::endl;

  Grid_finalize();
}
