    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_cg_prec.cc

    Copyright (C) 2015

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

using namespace std;
using namespace Grid;

#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#endif


NAMESPACE_BEGIN(Grid);
template<class Matrix,class Field>
  class SchurDiagMooeeOperatorParanoid :  public SchurOperatorBase<Field> {
 public:
    Matrix &_Mat;
    SchurDiagMooeeOperatorParanoid (Matrix &Mat): _Mat(Mat){};
    virtual  void Mpc      (const Field &in, Field &out) {
      Field tmp(in.Grid());
      tmp.Checkerboard() = !in.Checkerboard();
      //      std::cout <<" Mpc starting"<<std::endl;

      RealD nn = norm2(in); // std::cout <<" Mpc Prior to dslash norm is "<<nn<<std::endl;
      _Mat.Meooe(in,tmp);
      nn = norm2(tmp); //std::cout <<" Mpc Prior to Mooeinv "<<nn<<std::endl;
      _Mat.MooeeInv(tmp,out);
      nn = norm2(out); //std::cout <<" Mpc Prior to dslash norm is "<<nn<<std::endl;
      _Mat.Meooe(out,tmp);
      nn = norm2(tmp); //std::cout <<" Mpc Prior to Mooee "<<nn<<std::endl;
      _Mat.Mooee(in,out);
      nn = norm2(out); //std::cout <<" Mpc Prior to axpy "<<nn<<std::endl;
      axpy(out,-1.0,tmp,out);
    }
    virtual void MpcDag   (const Field &in, Field &out){
      Field tmp(in.Grid());
      //      std::cout <<" MpcDag starting"<<std::endl;
      RealD nn = norm2(in);// std::cout <<" MpcDag Prior to dslash norm is "<<nn<<std::endl;
      _Mat.MeooeDag(in,tmp);
      _Mat.MooeeInvDag(tmp,out);
      nn = norm2(out);// std::cout <<" MpcDag Prior to dslash norm is "<<nn<<std::endl;
      _Mat.MeooeDag(out,tmp);
      nn = norm2(tmp);// std::cout <<" MpcDag Prior to Mooee "<<nn<<std::endl;
      _Mat.MooeeDag(in,out);
      nn = norm2(out);// std::cout <<" MpcDag Prior to axpy "<<nn<<std::endl;
      axpy(out,-1.0,tmp,out);
    }
};

NAMESPACE_END(Grid);

int main (int argc, char ** argv)
{
  char hostname[HOST_NAME_MAX+1];
  gethostname(hostname, HOST_NAME_MAX+1);
  std::string host(hostname);
  
  Grid_init(&argc,&argv);

  const int Ls=12;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  GridCartesian         * UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         * FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian * FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);
  
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionD result(FGrid); result=Zero();
  LatticeGaugeFieldD Umu(UGrid);
  LatticeGaugeFieldF Umu_f(UGrid_f); 
  
  SU<Nc>::HotConfiguration(RNG4,Umu);

  precisionChange(Umu_f,Umu);
  
  RealD mass=0.1;
  RealD M5=1.8;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionF Ddwf_f(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);

  LatticeFermionD    src_o(FrbGrid);
  LatticeFermionD result_o(FrbGrid);
  LatticeFermionD result_o_2(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o.Checkerboard() = Odd;
  result_o = Zero();
  result_o_2.Checkerboard() = Odd;
  result_o_2 = Zero();

  SchurDiagMooeeOperatorParanoid<DomainWallFermionD,LatticeFermionD> HermOpEO(Ddwf);
  SchurDiagMooeeOperatorParanoid<DomainWallFermionF,LatticeFermionF> HermOpEO_f(Ddwf_f);

  int nsecs=600;
  if( GridCmdOptionExists(argv,argv+argc,"--seconds") ){
    std::string arg = GridCmdOptionPayload(argv,argv+argc,"--seconds");
    GridCmdOptionInt(arg,nsecs);
  }
  
  std::cout << GridLogMessage << "::::::::::::: Starting mixed CG for "<<nsecs <<" seconds" << std::endl;

  MixedPrecisionConjugateGradient<LatticeFermionD,LatticeFermionF> mCG(1.0e-8, 10000, 50, FrbGrid_f, HermOpEO_f, HermOpEO);
  double t1,t2,flops;
  double MdagMsiteflops = 1452; // Mobius (real coeffs)
  // CG overhead: 8 inner product, 4+8 axpy_norm, 4+4 linear comb (2 of)
  double CGsiteflops = (8+4+8+4+4)*Nc*Ns ;
  std:: cout << " MdagM site flops = "<< 4*MdagMsiteflops<<std::endl;
  std:: cout << " CG    site flops = "<< CGsiteflops <<std::endl;
  int iters;

  time_t now;
  time_t start = time(NULL);
  UGrid->Broadcast(0,(void *)&start,sizeof(start));

  FlightRecorder::ContinueOnFail = 0;
  FlightRecorder::PrintEntireLog = 0;
  FlightRecorder::ChecksumComms  = 1;
  FlightRecorder::ChecksumCommsSend=0;

  if(char *s=getenv("GRID_PRINT_ENTIRE_LOG"))  FlightRecorder::PrintEntireLog     = atoi(s);
  if(char *s=getenv("GRID_CHECKSUM_RECV_BUF")) FlightRecorder::ChecksumComms      = atoi(s);
  if(char *s=getenv("GRID_CHECKSUM_SEND_BUF")) FlightRecorder::ChecksumCommsSend  = atoi(s);

  int iter=0;
  do {
    if ( iter == 0 ) {
      FlightRecorder::SetLoggingMode(FlightRecorder::LoggingModeRecord);
    } else {
      FlightRecorder::SetLoggingMode(FlightRecorder::LoggingModeVerify);
    }
    std::cerr << "******************* SINGLE PRECISION SOLVE "<<iter<<std::endl;
    result_o = Zero();
    t1=usecond(); 
    mCG(src_o,result_o);
    t2=usecond(); 
    iters = mCG.TotalInnerIterations; //Number of inner CG iterations
    flops = MdagMsiteflops*4*FrbGrid->gSites()*iters;
    flops+= CGsiteflops*FrbGrid->gSites()*iters;
    std::cout << " SinglePrecision iterations/sec "<< iters/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " SinglePrecision GF/s "<< flops/(t2-t1)/1000.<<std::endl;
    std::cout << " SinglePrecision error count "<< FlightRecorder::ErrorCount()<<std::endl;

    assert(FlightRecorder::ErrorCount()==0);

    std::cout << " FlightRecorder is OK! "<<std::endl;
    iter ++;
    now = time(NULL); UGrid->Broadcast(0,(void *)&now,sizeof(now));
  } while (now < (start + nsecs/10) );
    
  std::cout << GridLogMessage << "::::::::::::: Starting double precision CG" << std::endl;
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  int i=0;
  do { 
    if ( i == 0 ) {
      FlightRecorder::SetLoggingMode(FlightRecorder::LoggingModeRecord);
    } else {
      FlightRecorder::SetLoggingMode(FlightRecorder::LoggingModeVerify);
    }
    std::cerr << "******************* DOUBLE PRECISION SOLVE "<<i<<std::endl;
    result_o_2 = Zero();
    t1=usecond(); 
    CG(HermOpEO,src_o,result_o_2);
    t2=usecond();
    iters = CG.IterationsToComplete;
    flops = MdagMsiteflops*4*FrbGrid->gSites()*iters; 
    flops+= CGsiteflops*FrbGrid->gSites()*iters;

    std::cout << " DoublePrecision iterations/sec "<< iters/(t2-t1)*1000.*1000.<<std::endl;
    std::cout << " DoublePrecision GF/s "<< flops/(t2-t1)/1000.<<std::endl;
    std::cout << " DoublePrecision error count "<< FlightRecorder::ErrorCount()<<std::endl;
    assert(FlightRecorder::ErrorCount()==0);
    std::cout << " FlightRecorder is OK! "<<std::endl;
    now = time(NULL); UGrid->Broadcast(0,(void *)&now,sizeof(now));
    i++;
  } while (now < (start + nsecs) );

  LatticeFermionD diff_o(FrbGrid);
  RealD diff = axpy_norm(diff_o, -1.0, result_o, result_o_2);

  std::cout << GridLogMessage << "::::::::::::: Diff between mixed and regular CG: " << diff << std::endl;
  assert(diff < 1e-4);
  
  Grid_finalize();
}
