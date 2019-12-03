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

#include <typeinfo>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// Very simple iterators for Eigen tensors
// The only way I could get these iterators to work is to put the begin() and end() functions in the Eigen namespace
// So if Eigen ever defines these, we'll have a conflict and have to change this
namespace Eigen {
  template <typename ET>
  inline typename std::enable_if<EigenIO::is_tensor<ET>::value, typename EigenIO::Traits<ET>::scalar_type *>::type
  begin( ET & et ) { return reinterpret_cast<typename Grid::EigenIO::Traits<ET>::scalar_type *>(et.data()); }
  template <typename ET>
  inline typename std::enable_if<EigenIO::is_tensor<ET>::value, typename EigenIO::Traits<ET>::scalar_type *>::type
  end( ET & et ) { return begin(et) + et.size() * EigenIO::Traits<ET>::count; }
}

/////////////////////////////////////////////////////////////
// Test creation of laplacian eigenvectors
/////////////////////////////////////////////////////////////

void test_Global(Application &application)
{
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start = 1100;
  globalPar.trajCounter.end   = 1120;
  globalPar.trajCounter.step  = 20;
  globalPar.runId             = "test";
  globalPar.graphFile         = "";
  globalPar.scheduleFile      = "";
  globalPar.saveSchedule      = "false";
  globalPar.parallelWriteMaxRetry      = -1;
  application.setPar(globalPar);
}

/////////////////////////////////////////////////////////////
// Create a random gauge with the correct name
/////////////////////////////////////////////////////////////

std::string test_Gauge(Application &application )
{
  std::string sGaugeName{ "gauge" };
  application.createModule<MGauge::Random>( sGaugeName );
  return sGaugeName;
}

/////////////////////////////////////////////////////////////
// Test creation of laplacian eigenvectors
/////////////////////////////////////////////////////////////

void test_LapEvec(Application &application)
{
  const char szModuleName[] = "LapEvec";
  test_Gauge( application );
  MDistil::LapEvecPar p;
  p.gauge = "gauge";
  p.Stout.steps = 3;
  p.Stout.rho = 0.2;
  p.Cheby.PolyOrder = 11;
  p.Cheby.alpha = 0.55;
  p.Cheby.beta = 35.5;
  p.Lanczos.Nvec = 5;
  p.Lanczos.Nk = 6;
  p.Lanczos.Np = 2;
  p.Lanczos.MaxIt = 1000;
  p.Lanczos.resid = 1e-2;
  p.Lanczos.IRLLog = 0;
  application.createModule<MDistil::LapEvec>(szModuleName,p);
}

/////////////////////////////////////////////////////////////
// Test creation Solver
/////////////////////////////////////////////////////////////

std::string SolverName( const char * pSuffix = nullptr ) {
  std::string sSolverName{ "CG" };
  if( pSuffix && pSuffix[0] ) {
    sSolverName.append( "_" );
    sSolverName.append( pSuffix );
  }
  return sSolverName;
}

std::string test_Solver(Application &application, const char * pSuffix = nullptr )
{
  std::string sActionName{ "DWF" };
  if( pSuffix && pSuffix[0] ) {
    sActionName.append( "_" );
    sActionName.append( pSuffix );
  }
  MAction::DWF::Par actionPar;
  actionPar.gauge = "gauge";
  actionPar.Ls    = 16;
  actionPar.M5    = 1.8;
  actionPar.mass  = 0.005;
  actionPar.boundary = "1 1 1 -1";
  actionPar.twist = "0. 0. 0. 0.";
  application.createModule<MAction::DWF>( sActionName, actionPar );
  MSolver::RBPrecCG::Par solverPar;
  solverPar.action       = sActionName;
  solverPar.residual     = 1.0e-2;
  solverPar.maxIteration = 10000;
  std::string sSolverName{ SolverName( pSuffix ) };
  application.createModule<MSolver::RBPrecCG>( sSolverName, solverPar );
  return sSolverName;
}

/////////////////////////////////////////////////////////////
// DistilParameters
/////////////////////////////////////////////////////////////

std::string test_DPar(Application &application) {
  MDistil::DistilParameters DPar;
  DPar.nvec = 5;
  DPar.nnoise = 1;
  DPar.tsrc = 0;
  DPar.LI = 5;
  DPar.TI = 8;
  DPar.SI = 4;
  std::string sDParName{"DPar_l"};
  application.createModule<MDistil::DistilPar>(sDParName,DPar);
  return sDParName;
}
/////////////////////////////////////////////////////////////
// Noises
/////////////////////////////////////////////////////////////

std::string test_Noises(Application &application, const std::string &sNoiseBaseName ) {
  MDistil::NoisesPar NoisePar;
  NoisePar.DistilParams = "DPar_l";
  NoisePar.NoiseFileName = "noise";
  std::string sNoiseName{"noise"};
  application.createModule<MDistil::Noises>(sNoiseName,NoisePar);
  return sNoiseName;
}

/////////////////////////////////////////////////////////////
// Perambulators
/////////////////////////////////////////////////////////////

std::string PerambulatorName( const char * pszSuffix = nullptr )
{
  std::string sPerambulatorName{ "Peramb" };
  if( pszSuffix && pszSuffix[0] )
    sPerambulatorName.append( pszSuffix );
  return sPerambulatorName;
}

void test_LoadPerambulators( Application &application, const char * pszSuffix = nullptr )
{
  std::string sModuleName{ "Peramb_load" };
  MIO::LoadPerambulator::Par PerambPar;
  PerambPar.PerambFileName = "Peramb";
  PerambPar.DistilParams = "DPar_l";
  test_Noises(application, sModuleName); // I want these written after solver stuff
  application.createModule<MIO::LoadPerambulator>( sModuleName, PerambPar );
}

void test_Perambulators( Application &application, const char * pszSuffix = nullptr )
{
  std::string sModuleName{ PerambulatorName( pszSuffix ) };
  // Perambulator parameters
  MDistil::Perambulator::Par PerambPar;
  PerambPar.lapevec = "LapEvec";
  PerambPar.PerambFileName = sModuleName;
  PerambPar.solver = test_Solver( application, pszSuffix );
  PerambPar.DistilParams = "DPar_l";
  PerambPar.noise = "noise";
  test_Noises(application, sModuleName); // I want these written after solver stuff
  application.createModule<MDistil::Perambulator>( sModuleName, PerambPar );
}

/////////////////////////////////////////////////////////////
// DistilVectors
/////////////////////////////////////////////////////////////

void test_DistilVectors(Application &application, const char * pszSuffix = nullptr, const char * pszNvec = nullptr )
{
  std::string sModuleName{"DistilVecs"};
  if( pszSuffix )
    sModuleName.append( pszSuffix );
  std::string sPerambName{"Peramb"};
  if( pszSuffix )
    sPerambName.append( pszSuffix );
  MDistil::DistilVectors::Par DistilVecPar;
  DistilVecPar.noise = "noise";
  DistilVecPar.rho = "rho";
  DistilVecPar.phi = "phi";
  DistilVecPar.perambulator = sPerambName;
  DistilVecPar.lapevec = "LapEvec";
  DistilVecPar.DistilParams = "DPar_l";
  application.createModule<MDistil::DistilVectors>(sModuleName,DistilVecPar);
}

/////////////////////////////////////////////////////////////
// MesonSink
/////////////////////////////////////////////////////////////

void test_MesonSink(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  //A2AMesonFieldPar.left="Peramb_unsmeared_sink";
  A2AMesonFieldPar.left="g5phi";
  A2AMesonFieldPar.right="Peramb_unsmeared_sink";
  A2AMesonFieldPar.output="DistilFields";
  A2AMesonFieldPar.gammas="Identity";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonSink",A2AMesonFieldPar);
}
/////////////////////////////////////////////////////////////
// MesonFields
/////////////////////////////////////////////////////////////

void test_MesonField(Application &application, const char * pszFileSuffix,
                     const char * pszObjectLeft = nullptr, const char * pszObjectRight = nullptr )
{
  // DistilVectors parameters
  if( pszObjectLeft == nullptr )
    pszObjectLeft = pszFileSuffix;
  if( pszObjectRight == nullptr )
    pszObjectRight = pszObjectLeft;
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="";
  A2AMesonFieldPar.right=A2AMesonFieldPar.left;
  A2AMesonFieldPar.left.append( pszObjectLeft );
  A2AMesonFieldPar.right.append( pszObjectRight );
  A2AMesonFieldPar.output="MesonSinks";
  A2AMesonFieldPar.output.append( pszFileSuffix );
  A2AMesonFieldPar.gammas="Identity";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  std::string sObjectName{"DistilMesonField"};
  sObjectName.append( pszFileSuffix );
  application.createModule<MContraction::A2AMesonField>(sObjectName, A2AMesonFieldPar);
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
  
  switch(iTestNum) {
    default: // 0
      LOG(Message) << "Computing Meson 2pt-function" << std::endl;
      test_Global( application );
      test_LapEvec( application );
      test_DPar( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_MesonField( application, "Phi", "phi" );
      test_MesonField( application, "Rho", "rho" );
      break;
    case 1:
      LOG(Message) << "Computing Meson 2pt-function by loading perambulators" << std::endl;
      test_Global( application );
      test_LapEvec( application );
      test_DPar( application );
      test_LoadPerambulators( application );
      test_DistilVectors( application, "_load" );
      test_MesonField( application, "Phi", "phi" );
      test_MesonField( application, "Rho", "rho" );
      break;
  }
  // execution
  static const char XmlFileName[] = "test_distil.xml";
  application.saveParameterFile( XmlFileName );

  const Grid::Coordinate &lat{GridDefaultLatt()};
  if( lat.size() == 4 && lat[0] == 4 && lat[1] == 4 && lat[2] == 4 && lat[3] == 8 )
    application.run();
  else
    LOG(Warning) << "The parameters in " << XmlFileName << " are designed to run on a laptop usid --grid 4.4.4.8" << std::endl;
  
  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
