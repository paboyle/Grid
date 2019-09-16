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
  application.setPar(globalPar);
}

/////////////////////////////////////////////////////////////
// Create a random gauge with the correct name
/////////////////////////////////////////////////////////////

std::string test_Gauge(Application &application, const char * pszBaseName )
{
  std::string sGaugeName{ pszBaseName };
  sGaugeName.append( "_gauge" );
  application.createModule<MGauge::Random>( sGaugeName );
  return sGaugeName;
}

/////////////////////////////////////////////////////////////
// Test creation of laplacian eigenvectors
/////////////////////////////////////////////////////////////

void test_LapEvec(Application &application)
{
  const char szModuleName[] = "LapEvec";
  test_Gauge( application, szModuleName );
  MDistil::LapEvecPar p;
  p.Stout.steps = 3;
  p.Stout.rho = 0.2;
  p.Cheby.PolyOrder = 11;
  p.Cheby.alpha = 0.55;
  p.Cheby.beta = 12.5;
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
  actionPar.gauge = "LapEvec_gauge";
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
// Noises
/////////////////////////////////////////////////////////////

std::string test_Noises(Application &application, const std::string &sNoiseBaseName ) {
  // DistilVectors parameters
  MDistil::NoisesPar NoisePar;
  NoisePar.nnoise = 1;
  NoisePar.nvec = 5;
  std::string sNoiseName{sNoiseBaseName + "_noise"};
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
  std::string sModuleName{ PerambulatorName( pszSuffix ) };
  MIO::LoadPerambulator::Par PerambPar;
  PerambPar.PerambFileName = sModuleName;
  PerambPar.Distil.tsrc = 0;
  PerambPar.Distil.nnoise = 1;
  PerambPar.nvec = 5;
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
  PerambPar.Distil.tsrc = 0;
  PerambPar.Distil.nnoise = 1;
  PerambPar.nvec = 5;
  test_Noises(application, sModuleName); // I want these written after solver stuff
  application.createModule<MDistil::Perambulator>( sModuleName, PerambPar );
}

/////////////////////////////////////////////////////////////
// DistilVectors
/////////////////////////////////////////////////////////////

#define TEST_DISTIL_VECTORS_COMMON \
std::string sModuleName{"DistilVecs"}; \
if( pszSuffix ) \
  sModuleName.append( pszSuffix ); \
std::string sPerambName{"Peramb"}; \
if( pszSuffix ) \
  sPerambName.append( pszSuffix ); \
MDistil::DistilVectors::Par DistilVecPar; \
DistilVecPar.noise = sPerambName + "_noise"; \
DistilVecPar.perambulator = sPerambName; \
DistilVecPar.lapevec = "LapEvec"; \
DistilVecPar.tsrc = 0; \
if( pszNvec ) \
  DistilVecPar.nvec = pszNvec

#define TEST_DISTIL_VECTORS_COMMON_END \
application.createModule<MDistil::DistilVectors>(sModuleName,DistilVecPar)

void test_DistilVectors(Application &application, const char * pszSuffix = nullptr, const char * pszNvec = nullptr )
{
  TEST_DISTIL_VECTORS_COMMON;
  TEST_DISTIL_VECTORS_COMMON_END;
}

void test_DistilVectorsSS(Application &application, const char * pszSink, const char * pszSource,
                          const char * pszSuffix = nullptr, const char * pszNvec = nullptr )
{
  TEST_DISTIL_VECTORS_COMMON;
  if( pszSink )
    DistilVecPar.sink = pszSink;
  if( pszSource )
    DistilVecPar.source = pszSource;
  TEST_DISTIL_VECTORS_COMMON_END;
}

/////////////////////////////////////////////////////////////
// Multiple Perambulators
/////////////////////////////////////////////////////////////

void test_MultiPerambulators(Application &application)
{
  test_Perambulators( application, "5" );
  MDistil::PerambFromSolve::Par SolvePar;
  SolvePar.eigenPack="LapEvec";
  SolvePar.PerambFileName="Peramb2";
  SolvePar.solve = "Peramb5_unsmeared_sink";
  SolvePar.Distil.nnoise = 1;
  SolvePar.Distil.LI=5;
  SolvePar.Distil.SI=4;
  SolvePar.Distil.TI=8;
  SolvePar.nvec=5;
  SolvePar.nvec_reduced=2;
  SolvePar.LI_reduced=2;
  application.createModule<MDistil::PerambFromSolve>("Peramb2",SolvePar);
  SolvePar.PerambFileName="Peramb3";
  SolvePar.nvec_reduced=3;
  SolvePar.LI_reduced=3;
  application.createModule<MDistil::PerambFromSolve>("Peramb3",SolvePar);

  test_DistilVectors( application, "2", "2" );
  test_DistilVectors( application, "3", "3" );
  test_DistilVectors( application, "5", "5" );

  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecs2_rho";
  A2AMesonFieldPar.right="DistilVecs2_rho";
  A2AMesonFieldPar.output="MesonSinksRho2";
  A2AMesonFieldPar.gammas="Identity";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldRho2",A2AMesonFieldPar);
  A2AMesonFieldPar.left="DistilVecs2_phi";
  A2AMesonFieldPar.right="DistilVecs2_phi";
  A2AMesonFieldPar.output="MesonSinksPhi2";
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldPhi2",A2AMesonFieldPar);
  A2AMesonFieldPar.left="DistilVecs3_rho";
  A2AMesonFieldPar.right="DistilVecs3_rho";
  A2AMesonFieldPar.output="MesonSinksRho3";
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldRho3",A2AMesonFieldPar);
  A2AMesonFieldPar.left="DistilVecs3_phi";
  A2AMesonFieldPar.right="DistilVecs3_phi";
  A2AMesonFieldPar.output="MesonSinksPhi3";
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldPhi3",A2AMesonFieldPar);
  A2AMesonFieldPar.left="DistilVecs5_rho";
  A2AMesonFieldPar.right="DistilVecs5_rho";
  A2AMesonFieldPar.output="MesonSinksRho5";
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldRho5",A2AMesonFieldPar);
  A2AMesonFieldPar.left="DistilVecs5_phi";
  A2AMesonFieldPar.right="DistilVecs5_phi";
  A2AMesonFieldPar.output="MesonSinksPhi5";
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldPhi5",A2AMesonFieldPar);
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
  A2AMesonFieldPar.left="DistilVecs";
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

/////////////////////////////////////////////////////////////
// g5*unsmeared
/////////////////////////////////////////////////////////////

#ifdef DISTIL_PRE_RELEASE
void test_g5_sinks(Application &application)
{
  // DistilVectors parameters
  MDistil::g5_multiply::Par g5_multiplyPar;
  g5_multiplyPar.input="Peramb_unsmeared_sink";
  g5_multiplyPar.nnoise = 1;
  g5_multiplyPar.LI=5;
  g5_multiplyPar.Ns=4;
  g5_multiplyPar.Nt_inv=1;
  application.createModule<MDistil::g5_multiply>("g5phi",g5_multiplyPar);
}

/////////////////////////////////////////////////////////////
// BaryonFields - phiphiphi - efficient
/////////////////////////////////////////////////////////////

void test_BaryonFieldPhi2(Application &application)
{
  // DistilVectors parameters
  MDistil::BC2::Par BC2Par;
  BC2Par.one="DistilVecs_phi";
  BC2Par.two="DistilVecs_phi";
  BC2Par.three="DistilVecs_phi";
  BC2Par.output="BaryonFieldPhi2";
  BC2Par.parity=1;
  BC2Par.mom={"0 0 0"};
  application.createModule<MDistil::BC2>("BaryonFieldPhi2",BC2Par);
}
/////////////////////////////////////////////////////////////
// BaryonFields - rhorhorho - efficient
/////////////////////////////////////////////////////////////

void test_BaryonFieldRho2(Application &application)
{
  // DistilVectors parameters
  MDistil::BC2::Par BC2Par;
  BC2Par.one="DistilVecs_rho";
  BC2Par.two="DistilVecs_rho";
  BC2Par.three="DistilVecs_rho";
  BC2Par.output="BaryonFieldRho2";
  BC2Par.parity=1;
  BC2Par.mom={"0 0 0"};
  application.createModule<MDistil::BC2>("BaryonFieldRho2",BC2Par);
}
/////////////////////////////////////////////////////////////
// BaryonFields - phiphiphi
/////////////////////////////////////////////////////////////

void test_BaryonFieldPhi(Application &application)
{
  // DistilVectors parameters
  MDistil::BContraction::Par BContractionPar;
  BContractionPar.one="DistilVecs_phi";
  BContractionPar.two="DistilVecs_phi";
  BContractionPar.three="DistilVecs_phi";
  BContractionPar.output="BaryonFieldPhi";
  BContractionPar.parity=1;
  BContractionPar.mom={"0 0 0"};
  application.createModule<MDistil::BContraction>("BaryonFieldPhi",BContractionPar);
}
/////////////////////////////////////////////////////////////
// BaryonFields - rhorhorho
/////////////////////////////////////////////////////////////

void test_BaryonFieldRho(Application &application)
{
  // DistilVectors parameters
  MDistil::BContraction::Par BContractionPar;
  BContractionPar.one="DistilVecs_rho";
  BContractionPar.two="DistilVecs_rho";
  BContractionPar.three="DistilVecs_rho";
  BContractionPar.output="BaryonFieldRho";
  BContractionPar.parity=1;
  BContractionPar.mom={"0 0 0"};
  application.createModule<MDistil::BContraction>("BaryonFieldRho",BContractionPar);
}
/////////////////////////////////////////////////////////////
// BaryonContraction
/////////////////////////////////////////////////////////////

void test_Baryon2pt(Application &application)
{
  // DistilVectors parameters
  MDistil::Baryon2pt::Par Baryon2ptPar;
  Baryon2ptPar.inputL="BaryonFieldPhi";
  Baryon2ptPar.inputR="BaryonFieldRho";
  Baryon2ptPar.quarksL="uud";
  Baryon2ptPar.quarksR="uud";
  Baryon2ptPar.output="C2_baryon";
  application.createModule<MDistil::Baryon2pt>("C2_b",Baryon2ptPar);
}
#endif

/////////////////////////////////////////////////////////////
// emField
/////////////////////////////////////////////////////////////
void test_em(Application &application)
{
  MGauge::StochEm::Par StochEmPar;
  StochEmPar.gauge=PhotonR::Gauge::feynman;
  StochEmPar.zmScheme=PhotonR::ZmScheme::qedL;
  application.createModule<MGauge::StochEm>("Em",StochEmPar);
}

/////////////////////////////////////////////////////////////
// MesonA2ASlash
/////////////////////////////////////////////////////////////

void test_Aslash(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AAslashField::Par A2AAslashFieldPar;
  A2AAslashFieldPar.left="g5phi";
  //A2AAslashFieldPar.right="DistilVecs_phi";
  A2AAslashFieldPar.right="Peramb_unsmeared_sink";
  A2AAslashFieldPar.output="unsmeared_Aslash";
  A2AAslashFieldPar.emField={"Em"};
  A2AAslashFieldPar.cacheBlock=2;
  A2AAslashFieldPar.block=4;
  application.createModule<MContraction::A2AAslashField>("Aslash_field",A2AAslashFieldPar);
}

/////////////////////////////////////////////////////////////
// MesonA2ASlashSequential
/////////////////////////////////////////////////////////////

void test_AslashSeq(Application &application)
{
  // DistilVectors parameters
  MSolver::A2AAslashVectors::Par A2AAslashVectorsPar;
  A2AAslashVectorsPar.vector="PerambS_unsmeared_sink";
  A2AAslashVectorsPar.emField="Em";
  A2AAslashVectorsPar.solver="CG_s";
  A2AAslashVectorsPar.output="AslashSeq";
  application.createModule<MSolver::A2AAslashVectors>("Aslash_seq",A2AAslashVectorsPar);
}
/////////////////////////////////////////////////////////////
// Aslash_perambulators
/////////////////////////////////////////////////////////////
void test_PerambulatorsSolve(Application &application)
{
  // Perambulator parameters
  MDistil::PerambFromSolve::Par PerambFromSolvePar;
  PerambFromSolvePar.eigenPack="LapEvec";
  PerambFromSolvePar.solve="Aslash_seq";
  PerambFromSolvePar.PerambFileName="perambAslashS.bin";
  PerambFromSolvePar.Distil.tsrc = 0;
  PerambFromSolvePar.Distil.nnoise = 1;
  PerambFromSolvePar.nvec=5;
  application.createModule<MDistil::PerambFromSolve>("PerambAslashS",PerambFromSolvePar);
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

typedef Grid::Hadrons::MDistil::NamedTensor<Complex,3,sizeof(Real)> MyTensor;

template<typename T>
typename std::enable_if<Grid::EigenIO::is_tensor<T>::value && !Grid::Hadrons::MDistil::is_named_tensor<T>::value>::type
DebugShowTensor(T &x, const char * n, std::string * pIndexNames=nullptr)
{
  const MyTensor::Index s{x.size()};
  std::cout << n << ".size() = " << s << std::endl;
  std::cout << n << ".NumDimensions = " << x.NumDimensions << " (TensorBase)" << std::endl;
  std::cout << n << ".NumIndices = " << x.NumIndices << std::endl;
  const auto d{x.dimensions()};
  //std::cout << n << ".dimensions().size() = " << d.size() << std::endl;
  std::cout << "Dimensions are ";
  for(auto i = 0; i < x.NumDimensions ; i++)
    std::cout << "[" << d[i] << "]";
  std::cout << std::endl;
  MyTensor::Index SizeCalculated{1};
  std::cout << "Dimensions again";
  for(int i=0 ; i < x.NumDimensions ; i++ ) {
    std::cout << " : [" << i;
    if( pIndexNames )
      std::cout << ", " << pIndexNames[i];
    std::cout << "]=" << x.dimension(i);
    SizeCalculated *= d[i];
  }
  std::cout << std::endl;
  std::cout << "SizeCalculated = " << SizeCalculated << std::endl;\
  assert( SizeCalculated == s );
  // Initialise
  assert( x.NumDimensions == 3 );
  for( int i = 0 ; i < d[0] ; i++ )
    for( int j = 0 ; j < d[1] ; j++ )
      for( int k = 0 ; k < d[2] ; k++ ) {
        x(i,j,k) = std::complex<double>(SizeCalculated, -SizeCalculated);
        SizeCalculated--;
      }
  // Show raw data
  std::cout << "Data follow : " << std::endl;
  typename T::Scalar * p = x.data();
  for( auto i = 0 ; i < s ; i++ ) {
    if( i ) std::cout << ", ";
    std::cout << n << ".data()[" << i << "]=" << * p++;
  }
  std::cout << std::endl;
}

template<typename T>
typename std::enable_if<Grid::Hadrons::MDistil::is_named_tensor<T>::value>::type
DebugShowTensor(T &x, const char * n)
{
  DebugShowTensor( x.tensor, n, &x.IndexNames[0] );
}

// Test whether typedef and underlying types are the same

void DebugTestTypeEqualities(void)
{
  Real    r1;
  RealD   r2;
  double  r3;
  const std::type_info &tr1{typeid(r1)};
  const std::type_info &tr2{typeid(r2)};
  const std::type_info &tr3{typeid(r3)};
  if( tr1 == tr2 && tr2 == tr3 )
    std::cout << "r1, r2 and r3 are the same type" << std::endl;
  else
    std::cout << "r1, r2 and r3 are different types" << std::endl;
  std::cout << "r1 is a " << tr1.name() << std::endl;
  std::cout << "r2 is a " << tr2.name() << std::endl;
  std::cout << "r3 is a " << tr3.name() << std::endl;
  
  // These are the same
  Complex             c1;
  std::complex<Real>  c2;
  const std::type_info &tc1{typeid(c1)};
  const std::type_info &tc2{typeid(c2)};
  const std::type_info &tc3{typeid(SpinVector::scalar_type)};
  if( tc1 == tc2 && tc2 == tc3)
    std::cout << "c1, c2 and SpinVector::scalar_type are the same type" << std::endl;
  else
    std::cout << "c1, c2 and SpinVector::scalar_type are different types" << std::endl;
  std::cout << "c1                      is a " << tc1.name() << std::endl;
  std::cout << "c2                      is a " << tc2.name() << std::endl;
  std::cout << "SpinVector::scalar_type is a " << tc3.name() << std::endl;
  
  // These are the same
  SpinVector                              s1;
  iSpinVector<Complex >                   s2;
  iScalar<iVector<iScalar<Complex>, Ns> > s3;
  const std::type_info &ts1{typeid(s1)};
  const std::type_info &ts2{typeid(s2)};
  const std::type_info &ts3{typeid(s3)};
  if( ts1 == ts2 && ts2 == ts3 )
    std::cout << "s1, s2 and s3 are the same type" << std::endl;
  else
    std::cout << "s1, s2 and s3 are different types" << std::endl;
  std::cout << "s1 is a " << ts1.name() << std::endl;
  std::cout << "s2 is a " << ts2.name() << std::endl;
  std::cout << "s3 is a " << ts3.name() << std::endl;
  
  // These are the same
  SpinColourVector            sc1;
  iSpinColourVector<Complex > sc2;
  const std::type_info &tsc1{typeid(sc1)};
  const std::type_info &tsc2{typeid(sc2)};
  if( tsc1 == tsc2 )
    std::cout << "sc1 and sc2 are the same type" << std::endl;
  else
    std::cout << "sc1 and sc2 are different types" << std::endl;
  std::cout << "sc1 is a " << tsc1.name() << std::endl;
  std::cout << "sc2 is a " << tsc2.name() << std::endl;
}

bool DebugEigenTest()
{
  {
    Eigen::TensorFixedSize<std::complex<double>,Eigen::Sizes<3,4,5>> x;
    DebugShowTensor(x, "fixed");
  }
  const char pszTestFileName[] = "test_tensor.bin";
  std::array<std::string,3> as={"Alpha", "Beta", "Gamma"};
  MyTensor x(as, 2,1,4);
  DebugShowTensor(x, "x");
  x.write(pszTestFileName);
  // Test initialisation of an array of strings
  for( auto a : as )
    std::cout << a << std::endl;
  Grid::Hadrons::MDistil::NamedTensor<Complex,3,sizeof(Real)> p{as,2,7,2};
  DebugShowTensor(p, "p");
  std::cout << "p.IndexNames follow" << std::endl;
  for( auto a : p.IndexNames )
    std::cout << a << std::endl;

  // Now see whether we can read a tensor back
  std::array<std::string,3> Names2={"Alpha", "Gamma", "Delta"};
  MyTensor y(Names2, 2,4,1);
  y.read(pszTestFileName);
  DebugShowTensor(y, "y");

  // Now see whether we can read a tensor back from an hdf5 file
  const char * pszFileName = "test";
  y.write(pszFileName);
  {
    MyTensor z;
    const char * pszName = "z1";
    DebugShowTensor(z, pszName);
    z.read(pszFileName);
    DebugShowTensor(z, pszName);
  }
  {
    MyTensor z(Names2,2,0,0);
    const char * pszName = "z2";
    DebugShowTensor(z, pszName);
    z.read(pszFileName);
    DebugShowTensor(z, pszName);
  }
  {
    // Now see whether we can read a tensor back from an xml file
    const char * pszXmlName = "test.xml";
    {
      XmlWriter w(pszXmlName);
      y.write<XmlWriter>(w);
    }
    MyTensor z;
    const char * pszName = "xml1";
    DebugShowTensor(z, pszName);
    XmlReader r(pszXmlName);
    z.read<XmlReader>(r);
    DebugShowTensor(z, pszName);
  }
  if((0)) // The following tests would fail
  {
    MyTensor z(Names2,2,0,78);
    //std::array<std::string,3> NamesBad={"Alpha", "Gamma", "Kilo"};
    //MyTensor z(NamesBad);
    const char * pszName = "zFail";
    DebugShowTensor(z, pszName);
    z.read(pszFileName);
    DebugShowTensor(z, pszName);
  }

  // Testing whether typedef produces the same type - yes it does

  DebugTestTypeEqualities();
  std::cout << std::endl;

  // How to access members of SpinColourVector
  SpinColourVector  sc;
  for( int s = 0 ; s < Ns ; s++ ) {
    auto cv{sc()(s)};
    iVector<Complex,Nc> c2{sc()(s)};
    std::cout << " cv is a " << typeid(cv).name() << std::endl;
    std::cout << " c2 is a " << typeid(c2).name() << std::endl;
    for( int c = 0 ; c < Nc ; c++ ) {
      Complex & z{cv(c)};
      std::cout << "  sc[spin=" << s << ", colour=" << c << "] = " << z << std::endl;
    }
  }
  // We could have removed the Lorentz index independently, but much easier to do as we do above
  iVector<iVector<Complex,Nc>,Ns> sc2{sc()};
  std::cout << "sc() is a " << typeid(sc()).name() << std::endl;
  std::cout << "sc2  is a " << typeid(sc2 ).name() << std::endl;

  // Or you can access elements directly
  std::complex<Real> z = sc()(0)(0);
  std::cout << "z = " << z << std::endl;
  sc()(3)(2) = std::complex<Real>{3.141,-3.141};
  std::cout << "sc()(3)(2) = " << sc()(3)(2) << std::endl;

  return true;
}

template <typename T>
void DebugGridTensorTest_print( int i )
{
  // std::cout << i << " : " << EigenIO::is_tensor<T>::value
  // << ", Rank " << EigenIO::Traits<T>::Rank
  // << ", count " << EigenIO::Traits<T>::count
  // << std::endl;
}

// begin() and end() are the minimum necessary to support range-for loops
// should really turn this into an iterator ...
template<typename T, int N>
class TestObject {
public:
  using value_type = T;
private:
  value_type * m_p;
public:
  TestObject() {
    m_p = reinterpret_cast<value_type *>(std::malloc(N * sizeof(value_type)));
  }
  ~TestObject() { std::free(m_p); }
  inline value_type * begin(void) { return m_p; }
  inline value_type * end(void) { return m_p + N; }
};

template<typename ET> typename std::enable_if<EigenIO::is_tensor<ET>::value>::type
dump_tensor(const ET & et, const char * psz = nullptr) {
  if( psz )
    std::cout << psz << ": ";
  else
    std::cout << "Unnamed tensor: ";
  Serializable::WriteMember( std::cout, et );
}

template <int Options>
void EigenSliceExample()
{
  std::cout << "Eigen example, Options = " << Options << std::endl;
  using T2 = Eigen::Tensor<int, 2, Options>;
  T2 a(4, 3);
  a.setValues({{0, 100, 200}, {300, 400, 500},
    {600, 700, 800}, {900, 1000, 1100}});
  std::cout << "a\n" << a << std::endl;
  dump_tensor( a, "a" );
  Eigen::array<typename T2::Index, 2> offsets = {0, 1};
  Eigen::array<typename T2::Index, 2> extents = {4, 2};
  T2 slice = a.slice(offsets, extents);
  std::cout << "slice\n" << slice << std::endl;
  dump_tensor( slice, "slice" );
  std::cout << "\n========================================" << std::endl;
}

template <int Options>
void EigenSliceExample2()
{
  using TestScalar = std::complex<float>;
  using T3 = Eigen::Tensor<TestScalar, 3, Options>;
  using T2 = Eigen::Tensor<TestScalar, 2, Options>;
  T3 a(2,3,4);

  std::cout << "Initialising a:";
  TestScalar f{ 0 };
  const TestScalar Inc{ 1, -1 };
  for( auto &c : a ) {
    c = f;
    f += Inc;
  }
  std::cout << std::endl;
  std::cout << "Validating   a (Eigen::" << ( ( Options & Eigen::RowMajor ) ? "Row" : "Col" ) << "Major):" << std::endl;
  f = 0;
  for( int i = 0 ; i < a.dimension(0) ; i++ )
    for( int j = 0 ; j < a.dimension(1) ; j++ )
      for( int k = 0 ; k < a.dimension(2) ; k++ ) {
        std::cout << " a(" << i << "," << j << "," << k << ")=" << a(i,j,k) << std::endl;
        assert( ( Options & Eigen::RowMajor ) == 0 || a(i,j,k) == f );
        f += Inc;
      }
  //std::cout << std::endl;
  //std::cout << "a initialised to:\n" << a << std::endl;
  dump_tensor( a, "a" );
  std::cout << std::endl;
  Eigen::array<typename T3::Index, 3> offsets = {0,1,1};
  Eigen::array<typename T3::Index, 3> extents = {1,2,2};
  T3 b;
  b = a.slice( offsets, extents );//.reshape(NewExtents);
  std::cout << "b = a.slice( offsets, extents ):\n" << b << std::endl;
  dump_tensor( b, "b" );
  T2 c(3,4);
  c = a.chip(0,1);
  std::cout << "c = a.chip(0,0):\n" << c << std::endl;
  dump_tensor( c, "c" );
  //T2 d = b.reshape(extents);
  //std::cout << "b.reshape(extents) is:\n" << d << std::endl;
  std::cout << "\n========================================" << std::endl;
}

void DebugFelixTensorTest( void )
{
  unsigned int Nmom = 2;
  unsigned int Nt = 2;
  unsigned int N_1 = 2;
  unsigned int N_2 = 2;
  unsigned int N_3 = 2;
  using BaryonTensorSet = Eigen::Tensor<Complex, 6, Eigen::RowMajor>;
  BaryonTensorSet BField3(Nmom,4,Nt,N_1,N_2,N_3);
  std::vector<Complex> Memory(Nmom * Nt * N_1 * N_2 * N_3 * 2);
  using BaryonTensorMap = Eigen::TensorMap<BaryonTensorSet>;
  BaryonTensorMap BField4 (&Memory[0], Nmom,4,Nt,N_1,N_2,N_3);

  EigenSliceExample<Eigen::RowMajor>();
  EigenSliceExample<0>();
  EigenSliceExample2<Eigen::RowMajor>();
  EigenSliceExample2<0>();
}

bool DebugGridTensorTest( void )
{
  DebugFelixTensorTest();
  typedef Complex t1;
  typedef iScalar<t1> t2;
  typedef iVector<t1, Ns> t3;
  typedef iMatrix<t1, Nc> t4;
  typedef iVector<iMatrix<t1,1>,4> t5;
  typedef iScalar<t5> t6;
  typedef iMatrix<t6, 3> t7;
  typedef iMatrix<iVector<iScalar<t7>,4>,2> t8;
  int i = 1;
  DebugGridTensorTest_print<t1>( i++ );
  DebugGridTensorTest_print<t2>( i++ );
  DebugGridTensorTest_print<t3>( i++ );
  DebugGridTensorTest_print<t4>( i++ );
  DebugGridTensorTest_print<t5>( i++ );
  DebugGridTensorTest_print<t6>( i++ );
  DebugGridTensorTest_print<t7>( i++ );
  DebugGridTensorTest_print<t8>( i++ );

  //using TOC7 = TestObject<std::complex<double>, 7>;
  using TOC7 = t7;
  TOC7 toc7;
  constexpr std::complex<double> Inc{1,-1};
  std::complex<double> Start{Inc};
  for( auto &x : toc7 ) {
    x = Start;
    Start += Inc;
  }
  i = 0;
  std::cout << "toc7:";
  for( auto x : toc7 ) std::cout << " [" << i++ << "]=" << x;
  std::cout << std::endl;

  //t2 o2;
  //auto a2 = TensorRemove(o2);
  //t3 o3;
  //t4 o4;
  //auto a3 = TensorRemove(o3);
  //auto a4 = TensorRemove(o4);
  
  return true;
}

bool ConvertPeramb(const char * pszSource, const char * pszDest) {
  Grid::Hadrons::MDistil::PerambTensor p(Hadrons::MDistil::PerambIndexNames);
  p.ReadBinary( pszSource );
  p.write(pszDest);
  return true;
}
#endif

int main(int argc, char *argv[])
{
#ifdef DEBUG
  // Debug only - test of Eigen::Tensor
  //if( DebugEigenTest() ) return 0;
  //if(DebugGridTensorTest()) return 0;
  //if(ConvertPeramb("PerambL_100_tsrc0.3000","PerambL_100_tsrc0.3000")) return 0;
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
    case 0:
      test_Global( application );
      test_LapEvec( application );
      break;
    case 1:
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      break;
    default: // 2
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      break;
    case 3:
      test_Global( application );
      test_LapEvec( application );
      test_LoadPerambulators( application );
      test_DistilVectors( application );
      break;
    case 4:
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_MesonField( application, "Phi", "_phi" );
      test_MesonField( application, "Rho", "_rho" );
      break;
    case 5:
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_Perambulators( application, "S" );
      test_DistilVectors( application, "S" );
      test_MesonField( application, "SPhi", "S_phi" );
      test_MesonField( application, "SRho", "S_rho" );
      break;
#ifdef DISTIL_PRE_RELEASE
    case 6: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_g5_sinks( application );
      test_MesonSink( application );
      break;
    case 7: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_BaryonFieldPhi( application );
      test_BaryonFieldRho( application );
      break;
#endif
    case 8: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_MesonField( application, "Phi", "_phi" );
      test_MesonField( application, "Rho", "_rho" );
      break;
#ifdef DISTIL_PRE_RELEASE
    case 9: // 3
      test_Global( application );
      test_Solver( application );
      test_Baryon2pt( application );
      break;
    case 10: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_g5_sinks( application );
      test_em( application );
      test_Aslash( application );
      break;
    case 11: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_BaryonFieldPhi2( application );
      test_BaryonFieldRho2( application );
      break;
#endif
    case 12: // 3
      test_Global( application );
      test_LapEvec( application );
      test_Perambulators( application, "S" );
      test_em( application );
      test_AslashSeq( application );
      test_PerambulatorsSolve( application );
      test_DistilVectorsSS( application, "AslashSeq", nullptr, "S" );
      test_MesonField( application, "AslashSeq" );
      break;
    case 13:
      test_Global( application );
      test_LapEvec( application );
      test_MultiPerambulators( application );
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
