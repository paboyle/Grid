/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_prec_change.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
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

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int Ls = 12;
  Coordinate latt4 = GridDefaultLatt();

  GridCartesian         * UGridD   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  GridCartesian         * FGridD   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridD);
  GridRedBlackCartesian * FrbGridD = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridD);

  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);

  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  
  std::cout << GridLogMessage << "Initialising 4d RNG" << std::endl;
  GridParallelRNG          RNG4(UGridD);  RNG4.SeedFixedIntegers(seeds4);
  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGridD);  RNG5.SeedFixedIntegers(seeds5);
  std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

  LatticeFermionD field_d(FGridD), tmp_d(FGridD);
  random(RNG5,field_d); tmp_d = field_d;

  LatticeFermionD2 field_d2(FGridF), tmp_d2(FGridF);
  precisionChange(field_d2, field_d); tmp_d2 = field_d2;

  LatticeFermionF field_f(FGridF), tmp_f(FGridF);
  precisionChange(field_f, field_d); tmp_f = field_f;

  int N = 500;

  double time_ds = 0, time_sd = 0;

  std::cout<<GridLogMessage << "Benchmarking single<->double original implementation (fields initially device-resident)" << std::endl;
  for(int i=0;i<N;i++){
    //We want to benchmark the typical scenario of both fields being device resident
    //To do this, invoke an operation that will open a device view and touch all sites
    //with a write operation that invalidates the CPU copy
    field_d = tmp_d;
    field_f = tmp_f;

    double start=usecond();
    precisionChangeOrig(field_d,field_f);
    double stop=usecond();
    time_sd += stop - start;

    field_d = tmp_d;
    field_f = tmp_f;

    start=usecond();
    precisionChangeOrig(field_f,field_d);
    stop=usecond();
    time_ds += stop - start;   
  }
  std::cout << "d->s " << time_ds/N << "us" << " s->d " << time_sd/N << "us" << std::endl;


  precisionChangeWorkspace wk_sp_to_dp(field_d.Grid(),field_f.Grid());
  precisionChangeWorkspace wk_dp_to_sp(field_f.Grid(),field_d.Grid());
  
  std::cout<<GridLogMessage << "Benchmarking single<->double with pregenerated workspace(fields initially device-resident)" << std::endl;
  time_sd = time_ds = 0;
  for(int i=0;i<N;i++){
    field_d = tmp_d;
    field_f = tmp_f;

    double start=usecond();
    precisionChange(field_d,field_f, wk_sp_to_dp);
    double stop=usecond();
    time_sd += stop - start;

    field_d = tmp_d;
    field_f = tmp_f;

    start=usecond();
    precisionChange(field_f,field_d, wk_dp_to_sp);
    stop=usecond();
    time_ds += stop - start;   
  }
  std::cout << "d->s " << time_ds/N << "us" << " s->d " << time_sd/N << "us" << std::endl;
  
  std::cout<<GridLogMessage << "Benchmarking single<->double with workspace generated on-the-fly (fields initially device-resident)" << std::endl;
  time_sd = time_ds = 0;
  for(int i=0;i<N;i++){
    field_d = tmp_d;
    field_f = tmp_f;

    double start=usecond();
    precisionChange(field_d,field_f);
    double stop=usecond();
    time_sd += stop - start;

    field_d = tmp_d;
    field_f = tmp_f;

    start=usecond();
    precisionChange(field_f,field_d);
    stop=usecond();
    time_ds += stop - start;

  }
  std::cout << "d->s " << time_ds/N << "us" << " s->d " << time_sd/N << "us" << std::endl;


  std::cout<<GridLogMessage << "Benchmarking single<->double2 (fields initially device-resident)" << std::endl;
  time_sd = time_ds = 0;
  for(int i=0;i<N;i++){
    field_d2 = tmp_d2;
    field_f = tmp_f;

    double start=usecond();
    precisionChangeFast(field_d2,field_f);
    double stop=usecond();
    time_sd += stop - start;

    field_d2 = tmp_d2;
    field_f = tmp_f;

    start=usecond();
    precisionChangeFast(field_f,field_d2);
    stop=usecond();
    time_ds += stop - start;
  }
  std::cout << "d->s " << time_ds/N << "us" << " s->d " << time_sd/N << "us" << std::endl;


  std::cout<<GridLogMessage << "Benchmarking single<->double2 through standard precisionChange call(fields initially device-resident) [NB: perf should be the same as the previous test!]" << std::endl;
  time_sd = time_ds = 0;
  for(int i=0;i<N;i++){
    field_d2 = tmp_d2;
    field_f = tmp_f;

    double start=usecond();
    precisionChange(field_d2,field_f);
    double stop=usecond();
    time_sd += stop - start;

    field_d2 = tmp_d2;
    field_f = tmp_f;

    start=usecond();
    precisionChange(field_f,field_d2);
    stop=usecond();
    time_ds += stop - start;
  }
  std::cout << "d->s " << time_ds/N << "us" << " s->d " << time_sd/N << "us" << std::endl;

  Grid_finalize();
}
