/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: ./tests/Test_gamma.cc

Copyright (C) 2015-2017

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Antonin Portelli <antonin.portelli@ed.ac.uk>

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
using namespace QCD;

static constexpr double                      tolerance = 1.0e-6;
static std::array<SpinMatrix, Gamma::nGamma> testAlgebra;

void print(const SpinMatrix &g)
{
  for(int i = 0; i < Ns; i++)
  {
    std::cout << GridLogMessage << "(";
    for(int j=0;j<Ns;j++){
      if ( abs(g()(i,j)()) == 0 ) {
        std::cout<< " 0";
      } else if ( abs(g()(i,j)() - Complex(0,1)) == 0){
        std::cout<< " i";
      } else if ( abs(g()(i,j)() + Complex(0,1)) == 0){
        std::cout<< "-i";
      } else if ( abs(g()(i,j)() - Complex(1,0)) == 0){
        std::cout<< " 1";
      } else if ( abs(g()(i,j)() + Complex(1,0)) == 0){
        std::cout<< "-1";
      }
      std::cout<<((j == Ns-1) ? ")" : "," );
    }
    std::cout << std::endl;
  }
  std::cout << GridLogMessage << std::endl;
}

void createTestAlgebra(void)
{
  std::array<SpinMatrix, 4> testg;
  SpinMatrix                testg5;
  const Complex             I(0., 1.), mI(0., -1.);
  
  testg[0] = zero;
  testg[0]()(0, 3) = I;
  testg[0]()(1, 2) = I;
  testg[0]()(2, 1) = mI;
  testg[0]()(3, 0) = mI;
  std::cout << GridLogMessage << "test GammaX= " << std::endl;
  print(testg[0]);
  testg[1] = zero;
  testg[1]()(0, 3) = -1.;
  testg[1]()(1, 2) = 1.;
  testg[1]()(2, 1) = 1.;
  testg[1]()(3, 0) = -1.;
  std::cout << GridLogMessage << "test GammaY= " << std::endl;
  print(testg[1]);
  testg[2] = zero;
  testg[2]()(0, 2) = I;
  testg[2]()(1, 3) = mI;
  testg[2]()(2, 0) = mI;
  testg[2]()(3, 1) = I;
  std::cout << GridLogMessage << "test GammaZ= " << std::endl;
  print(testg[2]);
  testg[3] = zero;
  testg[3]()(0, 2) = 1.;
  testg[3]()(1, 3) = 1.;
  testg[3]()(2, 0) = 1.;
  testg[3]()(3, 1) = 1.;
  std::cout << GridLogMessage << "test GammaT= " << std::endl;
  print(testg[3]);
  testg5 = testg[0]*testg[1]*testg[2]*testg[3];
  
#define DEFINE_TEST_G(g, exp)\
testAlgebra[Gamma::Algebra::g]        = exp;\
testAlgebra[Gamma::Algebra::Minus##g] = -exp;\
  
  DEFINE_TEST_G(Identity    , 1.);
  DEFINE_TEST_G(Gamma5      , testg5);
  DEFINE_TEST_G(GammaX      , testg[0]);
  DEFINE_TEST_G(GammaY      , testg[1]);
  DEFINE_TEST_G(GammaZ      , testg[2]);
  DEFINE_TEST_G(GammaT      , testg[3]);
  DEFINE_TEST_G(GammaXGamma5, testg[0]*testg5);
  DEFINE_TEST_G(GammaYGamma5, testg[1]*testg5);
  DEFINE_TEST_G(GammaZGamma5, testg[2]*testg5);
  DEFINE_TEST_G(GammaTGamma5, testg[3]*testg5);
  DEFINE_TEST_G(SigmaXY     , .5*(testg[0]*testg[1] - testg[1]*testg[0]));
  DEFINE_TEST_G(SigmaXZ     , .5*(testg[0]*testg[2] - testg[2]*testg[0]));
  DEFINE_TEST_G(SigmaXT     , .5*(testg[0]*testg[3] - testg[3]*testg[0]));
  DEFINE_TEST_G(SigmaYZ     , .5*(testg[1]*testg[2] - testg[2]*testg[1]));
  DEFINE_TEST_G(SigmaYT     , .5*(testg[1]*testg[3] - testg[3]*testg[1]));
  DEFINE_TEST_G(SigmaZT     , .5*(testg[2]*testg[3] - testg[3]*testg[2]));
  
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

void checkGamma(const Gamma::Algebra a, GridSerialRNG &rng)
{
  SpinVector v;
  SpinMatrix m, &testg = testAlgebra[a];
  Gamma      g(a);
  bool       pass = true;
  
  random(rng, v);
  random(rng, m);
  
  std::cout << GridLogMessage << "Checking " << Gamma::name[a] << ": ";
  std::cout << "vecmul ";
  test(g*v, testg*v);
  std::cout << "matlmul ";
  test(g*m, testg*m);
  std::cout << "matrmul ";
  test(m*g, m*testg);
  std::cout << std::endl;
}

void checkProd(const Gamma::Algebra a, const Gamma::Algebra b)
{
  SpinMatrix gm, testg = testAlgebra[a]*testAlgebra[b];
  Gamma      g = Gamma(a)*Gamma(b);
  bool       pass = true;
  
  std::cout << GridLogMessage << "Checking " << Gamma::name[a] << " * "
            << Gamma::name[b] << ": ";
  gm = 1.0;
  gm = g*gm;
  test(gm, testg);
  std::cout << "(= " << Gamma::name[g.g] << ")" << std::endl;
}

void checkAdj(const Gamma::Algebra a)
{
  SpinMatrix gm, testg = adj(testAlgebra[a]);
  Gamma      g(adj(Gamma(a)));
  bool       pass = true;
  
  std::cout << GridLogMessage << "Checking adj(" << Gamma::name[a] << "): ";
  gm = 1.0;
  gm = g*gm;
  test(gm, testg);
  std::cout << "(= " << Gamma::name[g.g] << ")" << std::endl;
}

void checkProject(GridSerialRNG &rng)
{
  SpinVector     rv, recon, full;
  HalfSpinVector hsp, hsm;
  
  random(rng, rv);
  
#define CHECK_PROJ(dir, gamma)\
std::cout << GridLogMessage << "Checking " << #dir << " projector: ";\
spProj##dir(hsm,rv);\
spRecon##dir(recon,hsm);\
test(recon, rv + Gamma(Gamma::Algebra::gamma)*rv);\
std::cout << std::endl;
  
  CHECK_PROJ(Xp, GammaX);
  CHECK_PROJ(Yp, GammaY);
  CHECK_PROJ(Zp, GammaZ);
  CHECK_PROJ(Tp, GammaT);
  CHECK_PROJ(5p, Gamma5);
  CHECK_PROJ(Xm, MinusGammaX);
  CHECK_PROJ(Ym, MinusGammaY);
  CHECK_PROJ(Zm, MinusGammaZ);
  CHECK_PROJ(Tm, MinusGammaT);
  CHECK_PROJ(5m, MinusGamma5);
  
#undef CHECK_PROJ
}

void checkGammaL(const Gamma::Algebra a, GridSerialRNG &rng)
{
  SpinVector v;
  SpinMatrix m, &testg = testAlgebra[a], pl;
  GammaL     gl(a);
  bool       pass = true;
  
  random(rng, v);
  random(rng, m);
  
  pl = testAlgebra[Gamma::Algebra::Identity]
       - testAlgebra[Gamma::Algebra::Gamma5];
  std::cout << GridLogMessage << "Checking left-projected " << Gamma::name[a] << ": ";
  std::cout << "vecmul ";
  test(gl*v, testg*pl*v);
  std::cout << "matlmul ";
  test(gl*m, testg*pl*m);
  std::cout << "matrmul ";
  test(m*gl, m*testg*pl);
  std::cout << std::endl;
}

int main(int argc, char *argv[])
{
  Grid_init(&argc,&argv);
  
  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  
  GridCartesian Grid(latt_size,simd_layout,mpi_layout);
  GridSerialRNG sRNG;
  
  sRNG.SeedRandomDevice();
  
  std::cout << GridLogMessage << "======== Test algebra" << std::endl;
  createTestAlgebra();
  std::cout << GridLogMessage << "======== Multiplication operators check" << std::endl;
  for (int i = 0; i < Gamma::nGamma; ++i)
  {
    checkGamma(i, sRNG);
  }
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "======== Algebra multiplication table check" << std::endl;
  for (int i = 0; i < Gamma::nGamma; ++i)
  for (int j = 0; j < Gamma::nGamma; ++j)
  {
    checkProd(i, j);
  }
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "======== Adjoints check" << std::endl;
  for (int i = 0; i < Gamma::nGamma; ++i)
  {
    checkAdj(i);
  }
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "======== Spin projectors check" << std::endl;
  checkProject(sRNG);
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "======== Gamma-left matrices check" << std::endl;
  for (int i = 0; i < Gamma::nGamma; ++i)
  {
    checkGammaL(i, sRNG);
  }
  
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
