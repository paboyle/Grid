    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_quenched_update.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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

#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/WilsonLoops.h>
#include <qcd/utils/SUn.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt({8,8,8,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(latt, 
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());
  
  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  ///////////////////////////////
  // Configuration of known size
  ///////////////////////////////
  LatticeGaugeField Umu(grid);
  Umu=1.0; // Cold start

  // RNG set up for test
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // SU3 colour operatoions
  LatticeColourMatrix link(grid);
  LatticeColourMatrix staple(grid);
  int mu=0;

  // Apply heatbath to the link
  RealD beta=6.0;

  int subsets[2] = { Even, Odd};
  LatticeInteger one(rbGrid);  one = 1; // fill with ones
  LatticeInteger mask(grid); 

  for(int sweep=0;sweep<1000;sweep++){

    RealD plaq = ColourWilsonLoops::avgPlaquette(Umu);

    std::cout<<GridLogMessage<<"sweep "<<sweep<<" PLAQUETTE "<<plaq<<std::endl;

    for( int cb=0;cb<2;cb++ ) {

      one.checkerboard=subsets[cb];
      mask= zero;
      setCheckerboard(mask,one);

      //      std::cout<<GridLogMessage<<mask<<std::endl;
      for(int mu=0;mu<Nd;mu++){
	
	// Get Link and Staple term in action; must contain Beta and 
	// any other coeffs
	ColourWilsonLoops::Staple(staple,Umu,mu);

	link = PeekIndex<LorentzIndex>(Umu,mu);

	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {

	  // update Even checkerboard
	  SU3::SubGroupHeatBath(sRNG,pRNG,beta,link,staple,subgroup,20,mask);

	}

	PokeIndex<LorentzIndex>(Umu,link,mu);
	
	//reunitarise link;
	ProjectOnGroup(Umu);

      }

    }
    
  }

  Grid_finalize();
}




