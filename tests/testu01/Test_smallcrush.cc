    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_smallcrush.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// Wrap Grid's parallel RNG for testU01

extern "C" { 
#include "TestU01.h"
#include "gdef.h"
#include "unif01.h"
#include "ucarry.h"
#include "bbattery.h"
}

std::vector<std::ranlux48>      EngineRanlux;
std::vector<std::mt19937>       EngineMT;
std::vector<sitmo::prng_engine> EngineSitmo;

std::uniform_int_distribution<uint32_t> uid;

uint32_t GetU01Ranlux(void) {
  return uid(EngineRanlux[0]);
};
uint32_t GetU01MT(void) {
  return uid(EngineMT[0]);
};
uint32_t GetU01Sitmo(void) {
  return uid(EngineSitmo[0]);
};

typedef Grid::GridRNGbase::RngEngine RngEngine;

struct TestRNG { 
public:
  static GridParallelRNG *pRNG;
  static GridSerialRNG *sRNG;
  static GridBase *_grid;
  static RngEngine Eng;
  static uint64_t site;
  static uint64_t gsites;
  static char *name;

  static void Init(GridParallelRNG *_pRNG,GridSerialRNG *_sRNG,GridBase *grid) {
    pRNG = _pRNG;
    sRNG = _sRNG;
    _grid= grid;
    gsites= grid->_gsites;
    std::cout << "Init: Global sites is " <<gsites<<std::endl;
    site = 0;
  }
  static uint32_t GetU01(void) { 
    //    std::uniform_int_distribution<uint32_t> uid;
    uint32_t ret_val;
#if 0
    ret_val = sRNG->_uid[0](sRNG->_generators[0]);
#else
    ret_val = pRNG->GlobalU01(site);
    site=(site+1)%gsites;
#endif
    //    std::cout << "site "<<site <<" "<<std::hex << ret_val <<std::dec<<std::endl; 
    return ret_val;

  }
};

GridParallelRNG *TestRNG::pRNG;
GridSerialRNG   *TestRNG::sRNG;
GridBase        *TestRNG::_grid;
RngEngine        TestRNG::Eng;
uint64_t         TestRNG::site;
uint64_t         TestRNG::gsites;

#ifdef RNG_SITMO
char * TestRNG::name = (char *)"Grid_Sitmo";
#endif
#ifdef RNG_RANLUX
char * TestRNG::name = (char *)"Grid_ranlux48";
#endif
#ifdef RNG_MT19937
char * TestRNG::name = (char *)"Grid_mt19937";
#endif

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
     
  GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

  std::vector<int> seeds({1,2,3,4});
  std::seed_seq seq(seeds.begin(),seeds.end());

  EngineRanlux.push_back(std::ranlux48(seq));
  EngineMT.push_back(std::mt19937(seq));
  EngineSitmo.push_back(sitmo::prng_engine(seq));

  std::cout << GridLogMessage<< "Initialising Grid RNGs "<<std::endl; 
  GridParallelRNG           pRNG(&Grid);   
  pRNG.SeedFixedIntegers(std::vector<int>({43,12,7019,9}));
  GridSerialRNG           sRNG;
  sRNG.SeedFixedIntegers(std::vector<int>({102,12,99,15}));
  std::cout << GridLogMessage<< "Initialised Grid RNGs "<<std::endl; 

  TestRNG::Init(&pRNG,&sRNG,&Grid);
  std::cout << GridLogMessage<< "Grid RNG's are  "<< std::string(TestRNG::name) <<std::endl; 

  unif01_Gen * gen;

  //  gen = ulcg_CreateLCG (2147483647, 397204094, 0, 12345);
  //  bbattery_SmallCrush (gen);

  //  gen =  ucarry_CreateRanlux (48, 0x12345);
  //  bbattery_SmallCrush (gen);

  /*
  std::cout << GridLogMessage<< "Testing Standalone Ranlux" <<std::endl; 
  gen = unif01_CreateExternGenBits ((char *)"GridRanlux",GetU01Ranlux);
  bbattery_SmallCrush (gen);
  unif01_DeleteExternGenBits(gen);
  std::cout << GridLogMessage<< "Testing Standalone Ranlux is complete" <<std::endl; 
  */
  /*
  std::cout << GridLogMessage<< "Testing Standalone Mersenne Twister" <<std::endl; 
  gen = unif01_CreateExternGenBits ((char *)"GridMT",GetU01MT);
  bbattery_SmallCrush (gen);
  unif01_DeleteExternGenBits(gen);
  std::cout << GridLogMessage<< "Testing Standalone Mersenne Twister is complete" <<std::endl; 

  std::cout << GridLogMessage<< "Testing Standalone Sitmo" <<std::endl; 
  gen = unif01_CreateExternGenBits ((char *)"GridSitmo",GetU01Sitmo);
  bbattery_SmallCrush (gen);
  unif01_DeleteExternGenBits(gen);
  std::cout << GridLogMessage<< "Testing Standalone Sitmo is complete" <<std::endl; 
  */

  //  gen = unif01_CreateExternGenBits ((char *)"xorshift", xorshift);
  std::cout << GridLogMessage<< "Testing Grid "<< std::string(TestRNG::name) <<std::endl; 
  gen = unif01_CreateExternGenBits(TestRNG::name,TestRNG::GetU01);
  bbattery_SmallCrush (gen);
  std::cout << GridLogMessage<< "Testing Grid "<< std::string(TestRNG::name)<<" is complete" <<std::endl; 

  Grid_finalize();
}

