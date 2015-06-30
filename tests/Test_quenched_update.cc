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
							GridDefaultSimd(Nd,vComplexF::Nsimd()),
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

    std::cout<<"sweep "<<sweep<<" PLAQUETTE "<<plaq<<std::endl;

    for( int cb=0;cb<2;cb++ ) {

      one.checkerboard=subsets[cb];
      mask= zero;
      setCheckerboard(mask,one);

      //      std::cout<<mask<<std::endl;
      for(int mu=0;mu<Nd;mu++){
	
	// Get Link and Staple term in action; must contain Beta and 
	// any other coeffs
	ColourWilsonLoops::Staple(staple,Umu,mu);
	staple = adj(staple);

	link = PeekIndex<LorentzIndex>(Umu,mu);

	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {

	  // update Even checkerboard
	  SU3::SubGroupHeatBath(sRNG,pRNG,beta,link,staple,subgroup,20,mask);

	}

	PokeIndex<LorentzIndex>(Umu,link,mu);
	//reunitarise link;
	
      }

    }
    
  }

  Grid_finalize();
}




