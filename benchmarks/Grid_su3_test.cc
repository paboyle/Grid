#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

void su3_test_mult(LatticeColourMatrix &z, LatticeColourMatrix &x,LatticeColourMatrix &y)
{
  z=x*y;
}
