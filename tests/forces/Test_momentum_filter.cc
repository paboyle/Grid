    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_force.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>

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

//Get the mu-directected links on the upper boundary and the bulk remainder
template<typename Field>
void getLinksBoundaryBulk(Field &bound, Field &bulk,  Field &from, const Coordinate &latt_size){
  bound = Zero(); bulk = Zero();
  for(int mu=0;mu<Nd;mu++){
    LatticeInteger mucoor(bound.Grid());
    LatticeCoordinate(mucoor, mu);

    bound = where( mucoor == (Integer)(latt_size[mu] - 1), from, bound );
    bulk = where( mucoor != (Integer)(latt_size[mu] - 1), from, bulk );
  }
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  std::vector<int> seeds({1,2,3,4});

  GridParallelRNG          pRNG(&Grid);
  GridSerialRNG            sRNG;
  pRNG.SeedFixedIntegers(seeds);
  sRNG.SeedFixedIntegers(seeds);

  typedef PeriodicGimplR Gimpl;
  typedef WilsonGaugeAction<Gimpl> GaugeAction;
  typedef NoHirep Representation; //fundamental
  typedef NoSmearing<Gimpl> Smearing;
  typedef MinimumNorm2<Gimpl, Smearing> Omelyan;
  typedef Gimpl::Field Field;
  typedef MomentumFilterApplyPhase<Field> Filter;
  Filter filter(&Grid);
  
  //Setup a filter that disables link update on links passing through the global lattice boundary
  typedef Filter::LatticeLorentzScalarType MaskType;
  typedef Filter::LorentzScalarType MaskSiteType;

  MaskSiteType zero, one;
  for(int mu=0;mu<Nd;mu++){
    zero(mu)()() = 0.;
    one(mu)()() = 1.;
  }
  MaskType zeroField(&Grid), oneField(&Grid);
  zeroField = zero;
  oneField = one;

  
  filter.phase = oneField; //make every site 1.0

  //Zero mu-directed links at upper boundary
  for(int mu=0;mu<Nd;mu++){
    LatticeInteger mucoor(&Grid);
    LatticeCoordinate(mucoor, mu);
    
    filter.phase = where( mucoor == (Integer)(latt_size[mu] - 1) , zeroField, filter.phase );
  }

  //Start with a random gauge field
  Field U(&Grid);
  SU<Nc>::HotConfiguration(pRNG,U);

  //Get the original links on the bulk and boundary for later use
  Field Ubnd_orig(&Grid), Ubulk_orig(&Grid);
  getLinksBoundaryBulk(Ubnd_orig, Ubulk_orig, U, latt_size);

  ActionSet<Field,Representation> actions(1);  
  double beta=6;
  GaugeAction gauge_action(beta);
  actions[0].push_back(&gauge_action);

  Smearing smear;
  IntegratorParameters params(1,1.); //1 MD step
  Omelyan integrator(&Grid, params, actions, smear);
  
  integrator.setMomentumFilter(filter);

  integrator.refresh(U, sRNG, pRNG); //doesn't actually change the gauge field

  //Check the momentum is zero on the boundary
  const auto &P = integrator.getMomentum();
  Field Pbnd(&Grid), Pbulk(&Grid);
  getLinksBoundaryBulk(Pbnd, Pbulk, const_cast<Field&>(P), latt_size);

  RealD Pbnd_nrm = norm2(Pbnd); //expect zero
  std::cout << GridLogMessage << "After refresh, norm2 of mu-directed conjugate momentum on boundary is: " << Pbnd_nrm << " (expect 0)" << std::endl;
  RealD Pbulk_nrm = norm2(Pbulk); //expect non-zero
  std::cout << GridLogMessage << "After refresh, norm2 of bulk conjugate momentum is: " << Pbulk_nrm << " (expect non-zero)" << std::endl;

  //Evolve the gauge field
  integrator.integrate(U);

  //Check momentum is still zero on boundary
  getLinksBoundaryBulk(Pbnd, Pbulk, const_cast<Field&>(P), latt_size);
  
  Pbnd_nrm = norm2(Pbnd); //expect zero
  std::cout << GridLogMessage << "After integrate, norm2 of mu-directed conjugate momentum on boundary is: " << Pbnd_nrm << " (expect 0)" << std::endl;
  Pbulk_nrm = norm2(Pbulk); //expect non-zero
  std::cout << GridLogMessage << "After integrate, norm2 of bulk conjugate momentum is: " << Pbulk_nrm << " (expect non-zero)" << std::endl;

  //Get the new bulk and bound links
  Field Ubnd_new(&Grid), Ubulk_new(&Grid);
  getLinksBoundaryBulk(Ubnd_new, Ubulk_new, U, latt_size);

  Field Ubnd_diff = Ubnd_new - Ubnd_orig;
  Field Ubulk_diff = Ubulk_new - Ubulk_orig;

  RealD Ubnd_change = norm2( Ubnd_diff );
  RealD Ubulk_change = norm2( Ubulk_diff );
  std::cout << GridLogMessage << "After integrate, norm2 of change in mu-directed boundary links is : " << Ubnd_change << " (expect 0)" << std::endl;
  std::cout << GridLogMessage << "After integrate, norm2 of change in bulk links is : " << Ubulk_change << " (expect non-zero)" << std::endl;

  Grid_finalize();
}
