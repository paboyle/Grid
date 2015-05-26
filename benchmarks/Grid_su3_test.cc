#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


void su3_test_mult_routine(LatticeColourMatrix &z, LatticeColourMatrix &x,LatticeColourMatrix &y)
{
  mult(z,x,y);
}
