#include <iostream>
using namespace std;

#include "Grid/Grid.h"


typedef Grid::vRealD vReal;

typedef vector<vReal, Grid::alignedAllocator<vReal>> vReal1;


const int n_simd_lanes = vReal::Nsimd();


inline long reduced  (
    const long number )
{
  return (number + n_simd_lanes - 1) / n_simd_lanes;
}


int main (void)
{

  cout << "Hello!" << endl;
  cout << "n SIMD lanes is " << n_simd_lanes << endl;


  const long N = 1000000;

  for (long n = 0; n < N; n++)
  {

  }


  return (0);

}
