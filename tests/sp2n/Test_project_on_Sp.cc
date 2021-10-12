#include <Grid/Grid.h>

using namespace Grid;

int main (int argc, char **argv)
{
     Grid_init(&argc,&argv);
    
    
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();

    GridCartesian             Grid(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian     RBGrid(&Grid);
    
    
    LatticeGaugeField Umu(&Grid);
    LatticeColourMatrixD U(&Grid);
    
    std::vector<int> pseeds({1,2,3,4,5});
    std::vector<int> sseeds({6,7,8,9,10});
    GridParallelRNG  pRNG(&Grid); pRNG.SeedFixedIntegers(pseeds);
    GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);
    
    SU<Nc>::HotConfiguration(pRNG,Umu);
    U = PeekIndex<LorentzIndex>(Umu,0);
    
    
    U = ProjectOnSpGroup(U);
    
    Grid_finalize();


}
