/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Tests/Hadrons/Test_hadrons_distil.cc
 
 Copyright (C) 2015-2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>

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

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

/////////////////////////////////////////////////////////////
// Test creation of laplacian eigenvectors
/////////////////////////////////////////////////////////////

void test_Global(Application &application)
{
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start = 1500;
  globalPar.trajCounter.end   = 1520;
  globalPar.trajCounter.step  = 20;
  globalPar.runId             = "test";
  application.setPar(globalPar);
}

/////////////////////////////////////////////////////////////
// Test creation of laplacian eigenvectors
/////////////////////////////////////////////////////////////

void test_LapEvec(Application &application)
{
  const char szGaugeName[] = "gauge";
  // gauge field
  application.createModule<MGauge::Random>(szGaugeName);
  // Now make an instance of the LapEvec object
  MDistil::LapEvecPar p;
  p.gauge = szGaugeName;
  //p.EigenPackName = "ePack";
  //p.Distil.TI = 8;
  //p.Distil.LI = 3;
  //p.Distil.Nnoise = 2;
  //p.Distil.tSrc = 0;
  p.Stout.steps = 3;
  p.Stout.parm = 0.2;
  p.Cheby.PolyOrder = 11;
  p.Cheby.alpha = 0.3;
  p.Cheby.beta = 12.5;
  p.Lanczos.Nvec = 5;
  p.Lanczos.Nk = 6;
  p.Lanczos.Np = 2;
  p.Lanczos.MaxIt = 1000;
  p.Lanczos.resid = 1e-2;
  application.createModule<MDistil::LapEvec>("LapEvec",p);
}

/////////////////////////////////////////////////////////////
// Perambulators
/////////////////////////////////////////////////////////////

void test_Perambulators(Application &application)
{
  // PerambLight parameters
  MDistil::PerambLight::Par PerambPar;
  PerambPar.eigenPack="LapEvec";
  PerambPar.PerambFileName="peramb.bin";
  PerambPar.Distil.tsrc = 0;
  PerambPar.Distil.nnoise = 1;
  PerambPar.Distil.LI=5;
  PerambPar.Distil.SI=4;
  PerambPar.Distil.TI=8;
  PerambPar.nvec=5;
  PerambPar.Distil.Ns=4;
  PerambPar.Distil.Nt=8;
  PerambPar.Distil.Nt_inv=1;
  PerambPar.Solver.mass=0.005;
  PerambPar.Solver.M5=1.8;
  PerambPar.Ls=16;
  PerambPar.Solver.CGPrecision=1e-8;
  PerambPar.Solver.MaxIterations=10000;
  application.createModule<MDistil::PerambLight>("Peramb",PerambPar);
}
/////////////////////////////////////////////////////////////
// DistilVectors
/////////////////////////////////////////////////////////////

void test_DistilVectors(Application &application)
{
  // DistilVectors parameters
  MDistil::DistilVectors::Par DistilVecPar;
  DistilVecPar.noise="Peramb_noise";
  DistilVecPar.perambulator="Peramb_perambulator_light";
  DistilVecPar.eigenPack="LapEvec";
  DistilVecPar.tsrc = 0;
  DistilVecPar.nnoise = 1;
  DistilVecPar.LI=5;
  DistilVecPar.SI=4;
  DistilVecPar.TI=8;
  DistilVecPar.nvec=5;
  DistilVecPar.Ns=4;
  DistilVecPar.Nt=8;
  DistilVecPar.Nt_inv=1;
  application.createModule<MDistil::DistilVectors>("DistilVecs",DistilVecPar);
}
void test_PerambulatorsS(Application &application)
{
  // PerambLight parameters
  MDistil::PerambLight::Par PerambPar;
  PerambPar.eigenPack="LapEvec";
  PerambPar.PerambFileName="perambS.bin";
  PerambPar.Distil.tsrc = 0;
  PerambPar.Distil.nnoise = 1;
  PerambPar.Distil.LI=3;
  PerambPar.Distil.SI=4;
  PerambPar.Distil.TI=8;
  PerambPar.nvec=3;
  PerambPar.Distil.Ns=4;
  PerambPar.Distil.Nt=8;
  PerambPar.Distil.Nt_inv=1;
  PerambPar.Solver.mass=0.04; //strange mass???
  PerambPar.Solver.M5=1.8;
  PerambPar.Ls=16;
  PerambPar.Solver.CGPrecision=1e-8;
  PerambPar.Solver.MaxIterations=10000;
  application.createModule<MDistil::PerambLight>("PerambS",PerambPar);
}
/////////////////////////////////////////////////////////////
// DistilVectors
/////////////////////////////////////////////////////////////

void test_DistilVectorsS(Application &application)
{
  // DistilVectors parameters
  MDistil::DistilVectors::Par DistilVecPar;
  DistilVecPar.noise="PerambS_noise";
  DistilVecPar.perambulator="PerambS_perambulator_light";
  DistilVecPar.eigenPack="LapEvec";
  DistilVecPar.tsrc = 0;
  DistilVecPar.nnoise = 1;
  DistilVecPar.LI=3;
  DistilVecPar.SI=4;
  DistilVecPar.TI=8;
  DistilVecPar.nvec=3;
  DistilVecPar.Ns=4;
  DistilVecPar.Nt=8;
  DistilVecPar.Nt_inv=1;
  application.createModule<MDistil::DistilVectors>("DistilVecsS",DistilVecPar);
}
/////////////////////////////////////////////////////////////
// MesonSink
/////////////////////////////////////////////////////////////

void test_MesonSink(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="Peramb_unsmeared_sink";
  A2AMesonFieldPar.right="Peramb_unsmeared_sink";
  A2AMesonFieldPar.output="DistilFields";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonSink",A2AMesonFieldPar);
}
/////////////////////////////////////////////////////////////
// MesonFields
/////////////////////////////////////////////////////////////

void test_MesonFieldSL(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecsS_phi";
  //A2AMesonFieldPar.right="DistilVecs_rho";
  A2AMesonFieldPar.right="DistilVecs_phi";
  A2AMesonFieldPar.output="DistilFieldsS";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldS",A2AMesonFieldPar);
}
/////////////////////////////////////////////////////////////
// MesonFields
/////////////////////////////////////////////////////////////

void test_MesonField(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecs_phi";
  //A2AMesonFieldPar.right="DistilVecs_rho";
  A2AMesonFieldPar.right="DistilVecs_phi";
  A2AMesonFieldPar.output="MesonSinks";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonField",A2AMesonFieldPar);
}

bool bNumber( int &ri, const char * & pstr, bool bGobbleWhiteSpace = true )
{
  if( bGobbleWhiteSpace )
    while( std::isspace(static_cast<unsigned char>(*pstr)) )
      pstr++;
  const char * p = pstr;
  bool bMinus = false;
  char c = * p++;
  if( c == '+' )
    c = * p++;
  else if( c == '-' ) {
    bMinus = true;
    c = * p++;
  }
  int n = c - '0';
  if( n < 0 || n > 9 )
    return false;
  while( * p >= '0' && * p <= '9' ) {
    n = n * 10 + ( * p ) - '0';
    p++;
  }
  if( bMinus )
    n *= -1;
  ri = n;
  pstr = p;
  return true;
}

#ifdef DEBUG

typedef Grid::Hadrons::MDistil::NamedTensor<Complex,Real,3> MyTensor;

void DebugShowTensor(MyTensor &x, const char * n)
{
  const MyTensor::Index s{x.size()};
  std::cout << n << ".size() = " << s << std::endl;
  std::cout << n << ".NumDimensions = " << x.NumDimensions << " (TensorBase)" << std::endl;
  std::cout << n << ".NumIndices = " << x.NumIndices << std::endl;
  const MyTensor::Dimensions & d{x.dimensions()};
  std::cout << n << ".dimensions().size() = " << d.size() << std::endl;
  std::cout << "Dimensions are ";
  for(auto i : d ) std::cout << "[" << i << "]";
  std::cout << std::endl;
  MyTensor::Index SizeCalculated{1};
  std::cout << "Dimensions again";
  for(int i=0 ; i < d.size() ; i++ ) {
    std::cout << " : [" << i << ", " << x.IndexNames[i] << "]=" << d[i];
    SizeCalculated *= d[i];
  }
  std::cout << std::endl;
  std::cout << "SizeCalculated = " << SizeCalculated << std::endl;\
  assert( SizeCalculated == s );
  // Initialise
  assert( d.size() == 3 );
  for( int i = 0 ; i < d[0] ; i++ )
  for( int j = 0 ; j < d[1] ; j++ )
  for( int k = 0 ; k < d[2] ; k++ ) {
    x(i,j,k) = std::complex<double>(SizeCalculated, -SizeCalculated);
    SizeCalculated--;
  }
  // Show raw data
  std::cout << "Data follow : " << std::endl;
  Complex * p = x.data();
  for( auto i = 0 ; i < s ; i++ ) {
    if( i ) std::cout << ", ";
    std::cout << n << ".data()[" << i << "]=" << * p++;
  }
  std::cout << std::endl;
}

bool DebugEigenTest()
{
  const char pszTestFileName[] = "test_tensor.bin";
  std::array<std::string,3> as={"Alpha", "Beta", "Gamma"};
  MyTensor x(as, 2,1,4);
  DebugShowTensor(x, "x");
  x.WriteBinary(pszTestFileName);
  DebugShowTensor(x, "x");
  // Test initialisation of an array of strings
  for( auto a : as )
    std::cout << a << std::endl;
  Grid::Hadrons::MDistil::Perambulator<Complex,Real,3> p{as,2,7,2};
  DebugShowTensor(p, "p");
  std::cout << "p.IndexNames follow" << std::endl;
  for( auto a : p.IndexNames )
    std::cout << a << std::endl;
  // Now see whether we can read a tensor back
  std::array<std::string,3> a2={"Alpha", "Gamma", "Delta"};
  MyTensor y(a2, 2,4,1);
  y.ReadBinary(pszTestFileName);
  DebugShowTensor(y, "y");
  return true;
}
#endif

int main(int argc, char *argv[])
{
#ifdef DEBUG
  // Debug only - test of Eigen::Tensor
  std::cout << "sizeof(std::streamsize) = " << sizeof(std::streamsize) << std::endl;
  std::cout << "sizeof(Eigen::Index) = " << sizeof(Eigen::Index) << std::endl;
  //if( DebugEigenTest() ) return 0;
#endif

  // Decode command-line parameters. 1st one is which test to run
  int iTestNum = -1;

  for(int i = 1 ; i < argc ; i++ ) {
    std::cout << "argv[" << i << "]=\"" << argv[i] << "\"" << std::endl;
    const char * p = argv[i];
    if( * p == '/' || * p == '-' ) {
      p++;
      char c = * p++;
      switch(toupper(c)) {
        case 'T':
          if( bNumber( iTestNum, p ) ) {
            std::cout << "Test " << iTestNum << " requested";
            if( * p )
              std::cout << " (ignoring trailer \"" << p << "\")";
            std::cout << std::endl;
          }
          else
            std::cout << "Invalid test \"" << &argv[i][2] << "\"" << std::endl;
          break;
        default:
          std::cout << "Ignoring switch \"" << &argv[i][1] << "\"" << std::endl;
          break;
      }
    }
  }

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  LOG(Message) << "Grid initialized" << std::endl;
  
  // run setup ///////////////////////////////////////////////////////////////
  Application              application;

  // For now perform free propagator test - replace this with distillation test(s)
  LOG(Message) << "====== Creating xml for test " << iTestNum << " ======" << std::endl;
  //const unsigned int  nt    = GridDefaultLatt()[Tp];
  
  switch(iTestNum) {
    case 1:
      test_Global( application );
      test_LapEvec( application );
      break;
    case 2:
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      break;
    case 3: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      break;
    default: // 4
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_MesonField( application );
      break;
    case 5: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_PerambulatorsS( application );
      test_DistilVectorsS( application );
      test_MesonFieldSL( application );
      break;
    case 6: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_MesonSink( application );
      break;
  }
  LOG(Message) << "====== XML creation for test " << iTestNum << " complete ======" << std::endl;

  // execution
  application.saveParameterFile("test_hadrons_distil.xml");
  application.run();
  
  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
