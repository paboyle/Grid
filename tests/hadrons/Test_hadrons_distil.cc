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
#include <Grid/util/EigenUtil.h>

using namespace Grid;
using namespace Hadrons;

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
// Test creation Solver
/////////////////////////////////////////////////////////////

void test_SolverS(Application &application)
{
    std::string boundary = "1 1 1 -1";


        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 16;
        actionPar.M5    = 1.8;
        actionPar.mass  = 0.005;
        actionPar.boundary = boundary;
        actionPar.twist = "0. 0. 0. 0.";
        application.createModule<MAction::DWF>("DWF_s", actionPar);


        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_s";
        solverPar.residual     = 1.0e-7;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_s", solverPar);
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
  //p.ConfigFileDir="/home/dp008/dp008/dc-rich6/Scripts/ConfigsDeflQED/";
  //p.ConfigFileName="ckpoint_lat.3000";
  p.ConfigFileDir="/home/dp008/dp008/paboyle/A2A/run/";
  p.ConfigFileName="ckpoint_lat.IEEE64BIG.1100";
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
  PerambPar.ConfigFileDir="/home/dp008/dp008/paboyle/A2A/run/";
  PerambPar.ConfigFileName="ckpoint_lat.IEEE64BIG.1100";
  PerambPar.UniqueIdentifier="full_dilution";
  PerambPar.solver="CG_s";
  PerambPar.Distil.tsrc = 0;
  PerambPar.Distil.nnoise = 1;
  PerambPar.Distil.LI=5;
  PerambPar.Distil.SI=4;
  PerambPar.Distil.TI=8;
  PerambPar.nvec=5;
  PerambPar.Distil.Ns=4;
  PerambPar.Distil.Nt=8;
  PerambPar.Distil.Nt_inv=1;
  //PerambPar.Solver.mass=0.005;
  //PerambPar.Solver.M5=1.8;
  //PerambPar.Ls=16;
  //PerambPar.Solver.CGPrecision=1e-8;
  //PerambPar.Solver.MaxIterations=10000;
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
  PerambPar.ConfigFileDir="/home/dp008/dp008/paboyle/A2A/run/";
  PerambPar.ConfigFileName="ckpoint_lat.IEEE64BIG.1100";
  PerambPar.UniqueIdentifier="full_dilution";
  PerambPar.solver="CG_s";
  PerambPar.Distil.tsrc = 0;
  PerambPar.Distil.nnoise = 1;
  PerambPar.Distil.LI=5;
  PerambPar.Distil.SI=4;
  PerambPar.Distil.TI=8;
  PerambPar.nvec=5;
  PerambPar.Distil.Ns=4;
  PerambPar.Distil.Nt=8;
  PerambPar.Distil.Nt_inv=1;
  //PerambPar.Solver.mass=0.005; //strange mass???
  //PerambPar.Solver.M5=1.8;
  //PerambPar.Ls=16;
  //PerambPar.Solver.CGPrecision=1e-8;
  //PerambPar.Solver.MaxIterations=10000;
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
  DistilVecPar.LI=5;
  DistilVecPar.SI=4;
  DistilVecPar.TI=32;
  DistilVecPar.nvec=5;
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
  //A2AMesonFieldPar.left="Peramb_unsmeared_sink";
  A2AMesonFieldPar.left="g5phi";
  A2AMesonFieldPar.right="Peramb_unsmeared_sink";
  A2AMesonFieldPar.output="DistilFields";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonSink",A2AMesonFieldPar);
}
/////////////////////////////////////////////////////////////
// g5*unsmeared
/////////////////////////////////////////////////////////////

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
// MesonFields - phiphi
/////////////////////////////////////////////////////////////

void test_MesonField(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecs_phi";
  //A2AMesonFieldPar.right="DistilVecs_rho";
  A2AMesonFieldPar.right="DistilVecs_phi";
  A2AMesonFieldPar.output="MesonSinksPhi";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonField",A2AMesonFieldPar);
}
/////////////////////////////////////////////////////////////
// MesonFields - rhorho
/////////////////////////////////////////////////////////////

void test_MesonFieldRho(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecs_rho";
  //A2AMesonFieldPar.right="DistilVecs_rho";
  A2AMesonFieldPar.right="DistilVecs_rho";
  A2AMesonFieldPar.output="MesonSinksRho";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldRho",A2AMesonFieldPar);
}
/////////////////////////////////////////////////////////////
// MesonFields - rhorhoall
/////////////////////////////////////////////////////////////

void test_MesonFieldRhoAll(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecs_rho_all_tsrc";
  A2AMesonFieldPar.right="DistilVecs_rho_all_tsrc";
  A2AMesonFieldPar.output="MesonSinksRhoAll";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldRhoAll",A2AMesonFieldPar);
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
  // PerambLight parameters
  MDistil::PerambFromSolve::Par PerambFromSolvePar;
  PerambFromSolvePar.eigenPack="LapEvec";
  PerambFromSolvePar.solve="Aslash_seq";
  PerambFromSolvePar.PerambFileName="perambAslashS.bin";
  PerambFromSolvePar.Distil.tsrc = 0;
  PerambFromSolvePar.Distil.nnoise = 1;
  PerambFromSolvePar.Distil.LI=5;
  PerambFromSolvePar.Distil.SI=4;
  PerambFromSolvePar.Distil.TI=8;
  PerambFromSolvePar.nvec=5;
  PerambFromSolvePar.Distil.Ns=4;
  PerambFromSolvePar.Distil.Nt=8;
  PerambFromSolvePar.Distil.Nt_inv=1;
  application.createModule<MDistil::PerambFromSolve>("PerambAslashS",PerambFromSolvePar);
}
/////////////////////////////////////////////////////////////
// DistilVectors
/////////////////////////////////////////////////////////////

void test_DistilVectorsAslashSeq(Application &application)
{
  // DistilVectors parameters
  MDistil::DistilSink::Par DistilSinkPar;
  DistilSinkPar.perambulator="PerambAslashS";
  DistilSinkPar.eigenPack="LapEvec";
  DistilSinkPar.tsrc = 0;
  DistilSinkPar.nnoise = 1;
  DistilSinkPar.LI=5;
  DistilSinkPar.SI=4;
  DistilSinkPar.TI=8;
  DistilSinkPar.nvec=5;
  DistilSinkPar.Ns=4;
  DistilSinkPar.Nt=8;
  DistilSinkPar.Nt_inv=1;
  application.createModule<MDistil::DistilSink>("DistilVecsAslashSeq",DistilSinkPar);
}
/////////////////////////////////////////////////////////////
// MesonFields - aslaaaash
/////////////////////////////////////////////////////////////

void test_MesonFieldAslashSeq(Application &application)
{
  // DistilVectors parameters
  MContraction::A2AMesonField::Par A2AMesonFieldPar;
  A2AMesonFieldPar.left="DistilVecsAslashSeq";
  //A2AMesonFieldPar.right="DistilVecs_rho";
  A2AMesonFieldPar.right="DistilVecsAslashSeq";
  A2AMesonFieldPar.output="MesonSinksAslashSeq";
  A2AMesonFieldPar.gammas="all";
  A2AMesonFieldPar.mom={"0 0 0"};
  A2AMesonFieldPar.cacheBlock=2;
  A2AMesonFieldPar.block=4;
  application.createModule<MContraction::A2AMesonField>("DistilMesonFieldAslashSeq",A2AMesonFieldPar);
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
void DebugShowTensor(T &x, const char * n)
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
    std::cout << " : [" << i << /*", " << x.IndexNames[i] << */"]=" << x.dimension(i);
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
  x.WriteBinary(pszTestFileName);
  DebugShowTensor(x, "x");
  // Test initialisation of an array of strings
  for( auto a : as )
    std::cout << a << std::endl;
  Grid::Hadrons::MDistil::Perambulator<Complex,3,sizeof(Real)> p{as,2,7,2};
  DebugShowTensor(p, "p");
  std::cout << "p.IndexNames follow" << std::endl;
  for( auto a : p.IndexNames )
    std::cout << a << std::endl;
  // Now see whether we can read a tensor back
  std::array<std::string,3> Names2={"Alpha", "Gamma", "Delta"};
  MyTensor y(Names2, 2,4,1);
  y.ReadBinary(pszTestFileName);
  DebugShowTensor(y, "y");

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
  SequentialInit( a );
  std::cout << std::endl;
  //std::cout << "Validating   a:";
  float z = 0;
  for( int i = 0 ; i < a.dimension(0) ; i++ )
    for( int j = 0 ; j < a.dimension(1) ; j++ )
      for( int k = 0 ; k < a.dimension(2) ; k++ ) {
        TestScalar w{z, -z};
        //std::cout << " a(" << i << "," << j << "," << k << ")=" << w;
        assert( a(i,j,k) == w );
        z++;
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

  t2 o2;
  auto a2 = TensorRemove(o2);
  //t3 o3;
  //t4 o4;
  //auto a3 = TensorRemove(o3);
  //auto a4 = TensorRemove(o4);
  
  return true;
}
#endif

int main(int argc, char *argv[])
{
#ifdef DEBUG
  // Debug only - test of Eigen::Tensor
  std::cout << "sizeof(int) = " << sizeof(int)
  << ", sizeof(long) = " << sizeof(long)
  << ", sizeof(size_t) = " << sizeof(size_t)
  << ", sizeof(std::size_t) = " << sizeof(std::size_t)
  << ", sizeof(std::streamsize) = " << sizeof(std::streamsize)
  << ", sizeof(Eigen::Index) = " << sizeof(Eigen::Index)
  << ", sizeof(hsize_t) = " << sizeof(hsize_t)
  << ", sizeof(unsigned long long) = " << sizeof(unsigned long long)
  << std::endl;
  if( DebugEigenTest() ) return 0;
  if(DebugGridTensorTest()) return 0;
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
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      break;
    case 3: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      break;
    default: // 4
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_MesonField( application );
      test_MesonFieldRhoAll( application );
      break;
    case 5: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_PerambulatorsS( application );
      test_DistilVectorsS( application );
      test_MesonFieldSL( application );
      break;
    case 6: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_g5_sinks( application );
      test_MesonSink( application );
      break;
    case 7: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_BaryonFieldPhi( application );
      test_BaryonFieldRho( application );
      break;
    case 8: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_MesonField( application );
      break;
    case 9: // 3
      test_Global( application );
      test_SolverS( application );
      test_Baryon2pt( application );
      break;
    case 10: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_g5_sinks( application );
      test_em( application );
      test_Aslash( application );
      break;
    case 11: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_Perambulators( application );
      test_DistilVectors( application );
      test_BaryonFieldPhi2( application );
      test_BaryonFieldRho2( application );
      break;
    case 12: // 3
      test_Global( application );
      test_SolverS( application );
      test_LapEvec( application );
      test_PerambulatorsS( application );
      test_em( application );
      test_AslashSeq( application );
      test_PerambulatorsSolve( application );
      test_DistilVectorsAslashSeq( application );
      test_MesonFieldAslashSeq( application );
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
