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
  PerambPar.tsrc = 0;
  PerambPar.nnoise = 1;
  PerambPar.LI=5;
  PerambPar.SI=4;
  PerambPar.TI=64;
  PerambPar.nvec=5;
  PerambPar.Ns=4;
  PerambPar.Nt=64;
  PerambPar.Nt_inv=1;
  PerambPar.mass=0.005;
  PerambPar.M5=1.8;
  PerambPar.Ls=16;
  PerambPar.CGPrecision=1e-8;
  PerambPar.MaxIterations=10000;
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
  DistilVecPar.LI=6;
  DistilVecPar.SI=4;
  DistilVecPar.TI=64;
  DistilVecPar.nvec=6;
  DistilVecPar.Ns=4;
  DistilVecPar.Nt=64;
  DistilVecPar.Nt_inv=1;
  application.createModule<MDistil::DistilVectors>("DistilVectorsInstance",DistilVecPar);
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

int main(int argc, char *argv[])
{
  // Decode command-line parameters. 1st one is which test to run
  int iTestNum = 2;

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
    default: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
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
