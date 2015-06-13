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

  std::vector<int> latt({4,4,4,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(latt, 
							GridDefaultSimd(Nd,vComplexF::Nsimd()),
							GridDefaultMpi());
  
  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(2)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  SU2::printGenerators();
  SU2::testGenerators();

  std::cout<<"*********************************************"<<std::endl;
  std::cout<<"* Generators for SU(3)"<<std::endl;
  std::cout<<"*********************************************"<<std::endl;
  SU3::printGenerators();
  SU3::testGenerators();

  //  std::cout<<"*********************************************"<<std::endl;
  //  std::cout<<"* Generators for SU(4)"<<std::endl;
  //  std::cout<<"*********************************************"<<std::endl;
  //  SU4::printGenerators();
  //  SU4::testGenerators();

  //  std::cout<<"*********************************************"<<std::endl;
  //  std::cout<<"* Generators for SU(5)"<<std::endl;
  //  std::cout<<"*********************************************"<<std::endl;
  //  SU5::printGenerators();
  //  SU5::testGenerators();

  ///////////////////////////////
  // Configuration of known size
  ///////////////////////////////
  NerscField header;
  std::string file("./ckpoint_lat.400");
  LatticeGaugeField Umu(grid);
  //  readNerscConfiguration(Umu,header,file);
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

  // Get Staple
  ColourWilsonLoops::Staple(staple,Umu,mu);
  // Get Link
  link = peekIndex<LorentzIndex>(Umu,mu);

  // Apply heatbath to the link
  RealD beta=6.0;
  int subgroup=0;
  int nhb=1;
  int trials=0;
  int fails=0;

  LatticeInteger one(rbGrid);  one = 1; // fill with ones
  LatticeInteger mask(grid);   mask= zero;
  one.checkerboard=Even;
  setCheckerboard(mask,one);

  // update Even checkerboard

  SU3::SubGroupHeatBath(sRNG,pRNG,beta,link,staple,subgroup,
			nhb,trials,fails,mask);



  Grid_finalize();
}


