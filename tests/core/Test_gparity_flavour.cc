/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_gparity_flavour.cc

Copyright (C) 2015-2017

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

using namespace Grid;

#ifdef ENABLE_GPARITY
static constexpr double                      tolerance = 1.0e-6;
static std::array<GparityFlavourMatrix, GparityFlavour::nSigma> testAlgebra;

void print(const GparityFlavourMatrix &g)
{
  for(int i = 0; i < Ngp; i++)
  {
    std::cout << GridLogMessage << "(";
    for(int j=0;j<Ngp;j++){
      if ( abs( g(i,j)()() ) == 0 ) {
        std::cout<< " 0";
      } else if ( abs(g(i,j)()() - Complex(0,1)) == 0){
        std::cout<< " i";
      } else if ( abs(g(i,j)()() + Complex(0,1)) == 0){
        std::cout<< "-i";
      } else if ( abs(g(i,j)()() - Complex(1,0)) == 0){
        std::cout<< " 1";
      } else if ( abs(g(i,j)()() + Complex(1,0)) == 0){
        std::cout<< "-1";
      }
      std::cout<<((j == Ngp-1) ? ")" : "," );
    }
    std::cout << std::endl;
  }
  std::cout << GridLogMessage << std::endl;
}

void createTestAlgebra(void)
{
  std::array<GparityFlavourMatrix, 3> testg;
  const Complex             I(0., 1.), mI(0., -1.);

  // 0 1
  // 1 0
  testg[0] = Zero();
  testg[0](0, 1)()() = 1.;
  testg[0](1, 0)()() = 1.;
  std::cout << GridLogMessage << "test SigmaX= " << std::endl;
  print(testg[0]);

  // 0 -i
  // i  0
  testg[1] = Zero();
  testg[1](0, 1)()() = mI;
  testg[1](1, 0)()() = I;
  std::cout << GridLogMessage << "test SigmaY= " << std::endl;
  print(testg[1]);

  // 1  0
  // 0 -1
  testg[2] = Zero();
  testg[2](0, 0)()() = 1.0;
  testg[2](1, 1)()() = -1.0;
  std::cout << GridLogMessage << "test SigmaZ= " << std::endl;
  print(testg[2]);

  
#define DEFINE_TEST_G(g, exp)\
testAlgebra[GparityFlavour::Algebra::g]        = exp; \
testAlgebra[GparityFlavour::Algebra::Minus##g] = -exp;
  
  DEFINE_TEST_G(SigmaX      , testg[0]);
  DEFINE_TEST_G(SigmaY      , testg[1]);
  DEFINE_TEST_G(SigmaZ      , testg[2]);
  DEFINE_TEST_G(Identity    , 1.);

  GparityFlavourMatrix pplus;
  pplus = 1.0;
  pplus = pplus + testg[1];
  pplus = pplus * 0.5;

  DEFINE_TEST_G(ProjPlus    , pplus);
  
  GparityFlavourMatrix pminus;
  pminus = 1.0;
  pminus = pminus - testg[1];
  pminus = pminus * 0.5;

  DEFINE_TEST_G(ProjMinus    , pminus);

#undef DEFINE_TEST_G
}

template <typename Expr>
void test(const Expr &a, const Expr &b)
{
  if (norm2(a - b) < tolerance)
  {
    std::cout << "[OK] ";
  }
  else
  {
    std::cout << "[fail]" << std::endl;
    std::cout << GridLogError << "a= " << a << std::endl;
    std::cout << GridLogError << "is different (tolerance= " << tolerance << ") from " << std::endl;
    std::cout << GridLogError << "b= " << b << std::endl;
    exit(EXIT_FAILURE);
  }
}

void checkSigma(const GparityFlavour::Algebra a, GridSerialRNG &rng)
{
  GparityFlavourVector v;
  GparityFlavourMatrix m, &testg = testAlgebra[a];
  GparityFlavour      g(a);
  
  random(rng, v);
  random(rng, m);
  
  std::cout << GridLogMessage << "Checking " << GparityFlavour::name[a] << ": ";
  std::cout << "vecmul ";
  test(g*v, testg*v);
  std::cout << "matlmul ";
  test(g*m, testg*m);
  std::cout << "matrmul ";
  test(m*g, m*testg);
  std::cout << std::endl;
}
#endif

int main(int argc, char *argv[])
{
  Grid_init(&argc,&argv);
#ifdef ENABLE_GPARITY  
  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  
  GridCartesian Grid(latt_size,simd_layout,mpi_layout);
  GridSerialRNG sRNG;
  
  sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  
  std::cout << GridLogMessage << "======== Test algebra" << std::endl;
  createTestAlgebra();
  std::cout << GridLogMessage << "======== Multiplication operators check" << std::endl;
  for (int i = 0; i < GparityFlavour::nSigma; ++i)
  {
    checkSigma(i, sRNG);
  }
  std::cout << GridLogMessage << std::endl;
#endif  
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
