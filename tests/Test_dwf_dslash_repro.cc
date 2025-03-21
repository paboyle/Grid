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

typedef LatticeFermionD FermionField;

int VerifyOnDevice(const FermionField &res, FermionField &ref)
{
  deviceVector<int> Fails(1);
  int * Fail = &Fails[0];
  int FailHost=0;
  
  typedef typename FermionField::vector_object vobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  const uint64_t NN = res.Grid()->oSites();

  acceleratorPut(*Fail,FailHost);

  accelerator_barrier();
  // Inject an error

  int injection=0;
  if(getenv("GRID_ERROR_INJECT")) injection=1;
  autoView(res_v,res,AcceleratorWrite);
  autoView(ref_v,ref,AcceleratorRead);
  if ( res.Grid()->ThisRank()== 0 )
  {
    if (((random()&0xF)==0)&&injection) {
      uint64_t sF = random()%(NN);
      int lane=0;
      printf("Error injection site %ld on rank %d\n",sF,res.Grid()->ThisRank());
      auto vv = acceleratorGet(res_v[sF]);
      double *dd = (double *)&vv;
      *dd=M_PI;
      acceleratorPut(res_v[sF],vv);
    }
  }

  accelerator_for( sF, NN, vobj::Nsimd(), {
#ifdef GRID_SIMT
      {
        int blane = acceleratorSIMTlane(vobj::Nsimd());
#else
      for(int blane;blane<vobj::Nsimd();blane++){
#endif
	vector_type *vtrr = (vector_type *)&res_v[sF];
	vector_type *vtrf = (vector_type *)&ref_v[sF];
	int words = sizeof(vobj)/sizeof(vector_type);
	
	for(int w=0;w<words;w++){
	  scalar_type rrtmp = getlane(vtrr[w], blane);
	  scalar_type rftmp = getlane(vtrf[w], blane);
	  if ( rrtmp != rftmp) {
	      *Fail=1;
	  }
	}
      }
  });

  FailHost = acceleratorGet(*Fail);

  return FailHost;
}
void PrintFails(const FermionField &res, FermionField &ref,uint64_t *ids)
{
  typedef typename FermionField::vector_object vobj;

  const int Nsimd=vobj::Nsimd();
  const uint64_t NN = res.Grid()->oSites();

  ///////////////////////////////
  // Pull back to host
  ///////////////////////////////
  autoView(res_v,res,CpuRead);
  autoView(ref_v,ref,CpuRead);
  
  std::vector<uint64_t> ids_host(NN*Nsimd);
  
  acceleratorCopyFromDevice(ids,&ids_host[0],NN*Nsimd*sizeof(uint64_t));

  //////////////////////////////////////////////////////////////
  // Redo check on host and print IDs
  //////////////////////////////////////////////////////////////
  
  for(int ss=0;ss< NN; ss++){				
      int sF = ss;
      for(int lane=0;lane<Nsimd;lane++){
	
	auto rr = extractLane(lane,res_v[sF]);
	auto rf = extractLane(lane,ref_v[sF]);
	uint64_t id = ids_host[lane+Nsimd*sF];
	//	std::cout << GridHostname()<<" id["<<sF<<"] lane "<<lane<<" id "<<id<<std::endl;
	for(int s=0;s<4;s++){
	  for(int c=0;c<3;c++){
	    if ( rr()(s)(c)!=rf()(s)(c) ) {
	      int subslice=(id>>0 )&0xFF;
	      int slice   =(id>>8 )&0xFF;
	      int eu      =(id>>16)&0xFF;
	      std::cout << GridHostname()<<" miscompare site "<<sF<<" "<<rr()(s)(c)<<" "<<rf()(s)(c)<<" EU "<<eu<<" slice "<<slice<<" subslice "<<subslice<<std::endl;
	    }
	  }
	}
      }
  };
  return;
}



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

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  LatticeFermionD    src(FGrid); random(RNG5,src);
  LatticeFermionD   junk(FGrid); random(RNG5,junk);

  LatticeFermionD result(FGrid); result=Zero();
  LatticeFermionD ref(FGrid); ref=Zero();
  
  SU<Nc>::HotConfiguration(RNG4,Umu);

  RealD mass=0.1;
  RealD M5=1.8;

  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  int nsecs=600;
  if( GridCmdOptionExists(argv,argv+argc,"--seconds") ){
    std::string arg = GridCmdOptionPayload(argv,argv+argc,"--seconds");
    GridCmdOptionInt(arg,nsecs);
  }
  
  std::cout << GridLogMessage << "::::::::::::: Job startup Barrier " << std::endl;
  UGrid->Barrier();
  std::cout << GridLogMessage << "::::::::::::: Job startup Barrier complete" << std::endl;

  std::cout << GridLogMessage << "::::::::::::: Starting DWF repro for "<<nsecs <<" seconds" << std::endl;

  time_t now;
  time_t start = time(NULL);
  UGrid->Broadcast(0,(void *)&start,sizeof(start));

  FlightRecorder::ContinueOnFail = 0;
  FlightRecorder::PrintEntireLog = 0;
  FlightRecorder::ChecksumComms  = 0;
  FlightRecorder::ChecksumCommsSend=0;

  if(char *s=getenv("GRID_PRINT_ENTIRE_LOG"))  FlightRecorder::PrintEntireLog     = atoi(s);
  if(char *s=getenv("GRID_CHECKSUM_RECV_BUF")) FlightRecorder::ChecksumComms      = atoi(s);
  if(char *s=getenv("GRID_CHECKSUM_SEND_BUF")) FlightRecorder::ChecksumCommsSend  = atoi(s);

  const uint64_t NN = FGrid->oSites()*vComplexD::Nsimd();
  
  deviceVector<uint64_t> ids_device(NN);
  uint64_t *ids = &ids_device[0];
  

  Ddwf.DhopComms(src,ref);
  Ddwf.DhopCalc(src,ref,ids);

  Ddwf.DhopComms(src,result);
  
  int iter=0;
  do {
    
    result=junk;

    Ddwf.DhopCalc(src,result,ids);

    if ( VerifyOnDevice(result, ref) ) {
      printf("Node %s Iter %d detected fails\n",GridHostname(),iter);
      PrintFails(result,ref,ids);
      //      std::cout << " Dslash "<<iter<<" is WRONG! "<<std::endl;
    }
    //else {
    //      printf("Node %s Iter %d detected NO fails\n",GridHostname(),iter);
    //      PrintFails(result,ref,ids);
    //      std::cout << " Dslash "<<iter<<" is OK! "<<std::endl;
    //}


    iter ++;
    now = time(NULL); UGrid->Broadcast(0,(void *)&now,sizeof(now));
  } while (now < (start + nsecs) );

  
  Grid_finalize();
}
