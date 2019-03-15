    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_memory_bandwidth.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
using namespace Grid::QCD;

typedef WilsonFermion5D<DomainWallVec5dImplR> WilsonFermion5DR;
typedef WilsonFermion5D<DomainWallVec5dImplF> WilsonFermion5DF;
typedef WilsonFermion5D<DomainWallVec5dImplD> WilsonFermion5DD;


std::vector<int> L_list;
std::vector<int> Ls_list;
std::vector<double> mflop_list;

double mflop_ref;
double mflop_ref_err;

int NN_global;

struct time_statistics{
  double mean;
  double err;
  double min;
  double max;

  void statistics(std::vector<double> v){
      double sum = std::accumulate(v.begin(), v.end(), 0.0);
      mean = sum / v.size();

      std::vector<double> diff(v.size());
      std::transform(v.begin(), v.end(), diff.begin(), [=](double x) { return x - mean; });
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      err = std::sqrt(sq_sum / (v.size()*(v.size() - 1)));

      auto result = std::minmax_element(v.begin(), v.end());
      min = *result.first;
      max = *result.second;
}
};

void comms_header(){
  std::cout <<GridLogMessage << " L  "<<"\t"<<" Ls  "<<"\t"
            <<std::setw(11)<<"bytes"<<"MB/s uni (err/min/max)"<<"\t\t"<<"MB/s bidi (err/min/max)"<<std::endl;
};

Gamma::Algebra Gmu [] = {
  Gamma::Algebra::GammaX,
  Gamma::Algebra::GammaY,
  Gamma::Algebra::GammaZ,
  Gamma::Algebra::GammaT
};
struct controls {
  int Opt;
  int CommsOverlap;
  Grid::CartesianCommunicator::CommunicatorPolicy_t CommsAsynch;
  //  int HugePages;
};

class Benchmark {
public:
  static void Decomposition (void ) {

    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Grid is setup to use "<<threads<<" threads"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage<<"Grid Default Decomposition patterns\n";
    std::cout<<GridLogMessage<<"\tOpenMP threads : "<<GridThread::GetThreads()<<std::endl;
    std::cout<<GridLogMessage<<"\tMPI tasks      : "<<GridCmdVectorIntToString(GridDefaultMpi())<<std::endl;
    std::cout<<GridLogMessage<<"\tvReal          : "<<sizeof(vReal )*8    <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vReal::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvRealF         : "<<sizeof(vRealF)*8    <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vRealF::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvRealD         : "<<sizeof(vRealD)*8    <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vRealD::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvComplex       : "<<sizeof(vComplex )*8 <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vComplex::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvComplexF      : "<<sizeof(vComplexF)*8 <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vComplexF::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvComplexD      : "<<sizeof(vComplexD)*8 <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vComplexD::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

  }

  static void Comms(void)
  {
    int Nloop=200;
    int nmu=0;
    int maxlat=32;

    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();

    for(int mu=0;mu<Nd;mu++) if (mpi_layout[mu]>1) nmu++;

    std::vector<double> t_time(Nloop);
    time_statistics timestat;

    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking threaded STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    comms_header();

    for(int lat=4;lat<=maxlat;lat+=4){
      for(int Ls=8;Ls<=8;Ls*=2){

	std::vector<int> latt_size  ({lat*mpi_layout[0],
	      lat*mpi_layout[1],
	      lat*mpi_layout[2],
	      lat*mpi_layout[3]});

	GridCartesian     Grid(latt_size,simd_layout,mpi_layout);
	RealD Nrank = Grid._Nprocessors;
	RealD Nnode = Grid.NodeCount();
	RealD ppn = Nrank/Nnode;

	std::vector<HalfSpinColourVectorD *> xbuf(8);
	std::vector<HalfSpinColourVectorD *> rbuf(8);
	Grid.ShmBufferFreeAll();
	for(int d=0;d<8;d++){
	  xbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
	  rbuf[d] = (HalfSpinColourVectorD *)Grid.ShmBufferMalloc(lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
	  bzero((void *)xbuf[d],lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
	  bzero((void *)rbuf[d],lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD));
	}

	int bytes=lat*lat*lat*Ls*sizeof(HalfSpinColourVectorD);
	int ncomm;
	double dbytes;
	std::vector<double> times(Nloop);
	for(int i=0;i<Nloop;i++){

	  double start=usecond();

	  dbytes=0;
	  ncomm=0;
#ifdef GRID_OMP
#pragma omp parallel for num_threads(Grid::CartesianCommunicator::nCommThreads)
#endif
	  for(int dir=0;dir<8;dir++){

	    double tbytes;
	    int mu =dir % 4;

	    if (mpi_layout[mu]>1 ) {
	        
	      int xmit_to_rank;
	      int recv_from_rank;
	      if ( dir == mu ) { 
		int comm_proc=1;
		Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      } else { 
		int comm_proc = mpi_layout[mu]-1;
		Grid.ShiftedRanks(mu,comm_proc,xmit_to_rank,recv_from_rank);
	      }
#ifdef GRID_OMP
	int tid = omp_get_thread_num(); 
#else 
        int tid = dir;
#endif
	      tbytes= Grid.StencilSendToRecvFrom((void *)&xbuf[dir][0], xmit_to_rank,
						 (void *)&rbuf[dir][0], recv_from_rank,
						 bytes,tid);
	  
#ifdef GRID_OMP
#pragma omp atomic
#endif
	      ncomm++;

#ifdef GRID_OMP
#pragma omp atomic
#endif
	      dbytes+=tbytes;
	    }
	  }
	  Grid.Barrier();
	  double stop=usecond();
	  t_time[i] = stop-start; // microseconds
	}

	timestat.statistics(t_time);
	//	for(int i=0;i<t_time.size();i++){
	//	  std::cout << i<<" "<<t_time[i]<<std::endl;
	//	}

	dbytes=dbytes*ppn;
	double xbytes    = dbytes*0.5;
	double rbytes    = dbytes*0.5;
	double bidibytes = dbytes;

	std::cout<<GridLogMessage << std::setw(4) << lat<<"\t"<<Ls<<"\t"
		 <<std::setw(11) << bytes<< std::fixed << std::setprecision(1) << std::setw(7)
		 <<std::right<< xbytes/timestat.mean<<"  "<< xbytes*timestat.err/(timestat.mean*timestat.mean)<< " "
		 <<xbytes/timestat.max <<" "<< xbytes/timestat.min  
		 << "\t\t"<<std::setw(7)<< bidibytes/timestat.mean<< "  " << bidibytes*timestat.err/(timestat.mean*timestat.mean) << " "
		 << bidibytes/timestat.max << " " << bidibytes/timestat.min << std::endl;

 
	
	    }
    }    

    return;
  }

  static void Memory(void)
  {
    const int Nvec=8;
    typedef Lattice< iVector< vReal,Nvec> > LatticeVec;
    typedef iVector<vReal,Nvec> Vec;

    std::vector<int> simd_layout = GridDefaultSimd(Nd,vReal::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();

    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking a*x + y bandwidth"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<"\t\t seconds"<< "\t\tGB/s / node"<<std::endl;
    std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;
  
    uint64_t NP;
    uint64_t NN;


  uint64_t lmax=48;
#define NLOOP (100*lmax*lmax*lmax*lmax/lat/lat/lat/lat)

    GridSerialRNG          sRNG;      sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    for(int lat=8;lat<=lmax;lat+=4){

      std::vector<int> latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol= latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      NP= Grid.RankCount();
      NN =Grid.NodeCount();

      Vec rn ; random(sRNG,rn);

      LatticeVec z(&Grid); z=rn;
      LatticeVec x(&Grid); x=rn;
      LatticeVec y(&Grid); y=rn;
      double a=2.0;

      uint64_t Nloop=NLOOP;

      double start=usecond();
      for(int i=0;i<Nloop;i++){
	z=a*x-y;
        x._odata[0]=z._odata[0]; // force serial dependency to prevent optimise away
        y._odata[4]=z._odata[4];
      }
      double stop=usecond();
      double time = (stop-start)/Nloop*1000;
     
      double flops=vol*Nvec*2;// mul,add
      double bytes=3.0*vol*Nvec*sizeof(Real);
      std::cout<<GridLogMessage<<std::setprecision(3) 
	       << lat<<"\t\t"<<bytes<<"   \t\t"<<bytes/time<<"\t\t"<<flops/time<<"\t\t"<<(stop-start)/1000./1000.
	       << "\t\t"<< bytes/time/NN <<std::endl;

    }
  };

  static double DWF5(int Ls,int L)
  {
    RealD mass=0.1;
    RealD M5  =1.8;

    double mflops;
    double mflops_best = 0;
    double mflops_worst= 0;
    std::vector<double> mflops_all;

    ///////////////////////////////////////////////////////
    // Set/Get the layout & grid size
    ///////////////////////////////////////////////////////
    int threads = GridThread::GetThreads();
    std::vector<int> mpi = GridDefaultMpi(); assert(mpi.size()==4);
    std::vector<int> local({L,L,L,L});

    GridCartesian         * TmpGrid   = SpaceTimeGrid::makeFourDimGrid(std::vector<int>({64,64,64,64}), 
								       GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
    NN_global=NN;
    uint64_t SHM=NP/NN;

    std::vector<int> internal;
    if      ( SHM == 1 )   internal = std::vector<int>({1,1,1,1});
    else if ( SHM == 2 )   internal = std::vector<int>({2,1,1,1});
    else if ( SHM == 4 )   internal = std::vector<int>({2,2,1,1});
    else if ( SHM == 8 )   internal = std::vector<int>({2,2,2,1});
    else assert(0);

    std::vector<int> nodes({mpi[0]/internal[0],mpi[1]/internal[1],mpi[2]/internal[2],mpi[3]/internal[3]});
    std::vector<int> latt4({local[0]*nodes[0],local[1]*nodes[1],local[2]*nodes[2],local[3]*nodes[3]});

    ///////// Welcome message ////////////
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "Benchmark DWF Ls vec on "<<L<<"^4 local volume "<<std::endl;
    std::cout<<GridLogMessage << "* Global volume  : "<<GridCmdVectorIntToString(latt4)<<std::endl;
    std::cout<<GridLogMessage << "* Ls             : "<<Ls<<std::endl;
    std::cout<<GridLogMessage << "* MPI ranks      : "<<GridCmdVectorIntToString(mpi)<<std::endl;
    std::cout<<GridLogMessage << "* Intranode      : "<<GridCmdVectorIntToString(internal)<<std::endl;
    std::cout<<GridLogMessage << "* nodes          : "<<GridCmdVectorIntToString(nodes)<<std::endl;
    std::cout<<GridLogMessage << "* Using "<<threads<<" threads"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    ///////// Lattice Init ////////////
    GridCartesian         * UGrid    = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    GridRedBlackCartesian * UrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    GridCartesian         * sUGrid   = SpaceTimeGrid::makeFourDimDWFGrid(latt4,GridDefaultMpi());
    GridRedBlackCartesian * sUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(sUGrid);
    GridCartesian         * sFGrid   = SpaceTimeGrid::makeFiveDimDWFGrid(Ls,UGrid);
    GridRedBlackCartesian * sFrbGrid = SpaceTimeGrid::makeFiveDimDWFRedBlackGrid(Ls,UGrid);

    ///////// RNG Init ////////////
    std::vector<int> seeds4({1,2,3,4});
    std::vector<int> seeds5({5,6,7,8});
    GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
    GridParallelRNG          RNG5(sFGrid);  RNG5.SeedFixedIntegers(seeds5);
    std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

    ///////// Source preparation ////////////
    LatticeFermion src   (sFGrid); random(RNG5,src);
    LatticeFermion tmp   (sFGrid);

    RealD N2 = 1.0/::sqrt(norm2(src));
    src = src*N2;
    
    LatticeGaugeField Umu(UGrid);  SU3::HotConfiguration(RNG4,Umu); 

    WilsonFermion5DR sDw(Umu,*sFGrid,*sFrbGrid,*sUGrid,*sUrbGrid,M5);
    LatticeFermion src_e (sFrbGrid);
    LatticeFermion src_o (sFrbGrid);
    LatticeFermion r_e   (sFrbGrid);
    LatticeFermion r_o   (sFrbGrid);
    LatticeFermion r_eo  (sFGrid);
    LatticeFermion err   (sFGrid);
    {

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd,src_o,src);

#if defined(AVX512) 
      const int num_cases = 6;
      std::string fmt("A/S ; A/O ; U/S ; U/O ; G/S ; G/O ");
#else
      const int num_cases = 4;
      std::string fmt("U/S ; U/O ; G/S ; G/O ");
#endif
      controls Cases [] = {
#ifdef AVX512
	{ QCD::WilsonKernelsStatic::OptInlineAsm , QCD::WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptInlineAsm , QCD::WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  },
#endif
	{ QCD::WilsonKernelsStatic::OptHandUnroll, QCD::WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptHandUnroll, QCD::WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptGeneric   , QCD::WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptGeneric   , QCD::WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  }
      }; 

      for(int c=0;c<num_cases;c++) {

	QCD::WilsonKernelsStatic::Comms = Cases[c].CommsOverlap;
	QCD::WilsonKernelsStatic::Opt   = Cases[c].Opt;
	CartesianCommunicator::SetCommunicatorPolicy(Cases[c].CommsAsynch);

	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
	if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
	if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
	if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

	int nwarm = 100;
	uint64_t ncall = 1000;

	double t0=usecond();
	sFGrid->Barrier();
	for(int i=0;i<nwarm;i++){
	  sDw.DhopEO(src_o,r_e,DaggerNo);
	}
	sFGrid->Barrier();
	double t1=usecond();

	sDw.ZeroCounters();
	time_statistics timestat;
	std::vector<double> t_time(ncall);
	for(uint64_t i=0;i<ncall;i++){
	  t0=usecond();
	  sDw.DhopEO(src_o,r_e,DaggerNo);
	  t1=usecond();
	  t_time[i] = t1-t0;
	}
	sFGrid->Barrier();
	
	double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
	double flops=(1344.0*volume)/2;
	double mf_hi, mf_lo, mf_err;

	timestat.statistics(t_time);
	mf_hi = flops/timestat.min;
	mf_lo = flops/timestat.max;
	mf_err= flops/timestat.min * timestat.err/timestat.mean;

	mflops = flops/timestat.mean;
	mflops_all.push_back(mflops);
	if ( mflops_best == 0   ) mflops_best = mflops;
	if ( mflops_worst== 0   ) mflops_worst= mflops;
	if ( mflops>mflops_best ) mflops_best = mflops;
	if ( mflops<mflops_worst) mflops_worst= mflops;

	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"sDeo mflop/s =   "<< mflops << " ("<<mf_err<<") " << mf_lo<<"-"<<mf_hi <<std::endl;
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"sDeo mflop/s per rank   "<< mflops/NP<<std::endl;
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"sDeo mflop/s per node   "<< mflops/NN<<std::endl;

	sDw.Report();

      }
      double robust = mflops_worst/mflops_best;;
      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
      std::cout<<GridLogMessage << L<<"^4 x "<<Ls<< " sDeo Best  mflop/s        =   "<< mflops_best << " ; " << mflops_best/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage << L<<"^4 x "<<Ls<< " sDeo Worst mflop/s        =   "<< mflops_worst<< " ; " << mflops_worst/NN<<" per node " <<std::endl;

      std::cout<<GridLogMessage <<std::setprecision(3)<< L<<"^4 x "<<Ls<< " Performance Robustness   =   "<< robust <<std::endl;
      std::cout<<GridLogMessage <<fmt << std::endl;
      std::cout<<GridLogMessage;

      for(int i=0;i<mflops_all.size();i++){
	std::cout<<mflops_all[i]/NN<<" ; " ;
      }
      std::cout<<std::endl;
      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    }
    return mflops_best;
  }

  static double DWF(int Ls,int L, double & robust)
  {
    RealD mass=0.1;
    RealD M5  =1.8;

    double mflops;
    double mflops_best = 0;
    double mflops_worst= 0;
    std::vector<double> mflops_all;

    ///////////////////////////////////////////////////////
    // Set/Get the layout & grid size
    ///////////////////////////////////////////////////////
    int threads = GridThread::GetThreads();
    std::vector<int> mpi = GridDefaultMpi(); assert(mpi.size()==4);
    std::vector<int> local({L,L,L,L});

    GridCartesian         * TmpGrid   = SpaceTimeGrid::makeFourDimGrid(std::vector<int>({64,64,64,64}), 
								       GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
    NN_global=NN;
    uint64_t SHM=NP/NN;

    std::vector<int> internal;
    if      ( SHM == 1 )   internal = std::vector<int>({1,1,1,1});
    else if ( SHM == 2 )   internal = std::vector<int>({2,1,1,1});
    else if ( SHM == 4 )   internal = std::vector<int>({2,2,1,1});
    else if ( SHM == 8 )   internal = std::vector<int>({2,2,2,1});
    else assert(0);

    std::vector<int> nodes({mpi[0]/internal[0],mpi[1]/internal[1],mpi[2]/internal[2],mpi[3]/internal[3]});
    std::vector<int> latt4({local[0]*nodes[0],local[1]*nodes[1],local[2]*nodes[2],local[3]*nodes[3]});

    ///////// Welcome message ////////////
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "Benchmark DWF on "<<L<<"^4 local volume "<<std::endl;
    std::cout<<GridLogMessage << "* Global volume  : "<<GridCmdVectorIntToString(latt4)<<std::endl;
    std::cout<<GridLogMessage << "* Ls             : "<<Ls<<std::endl;
    std::cout<<GridLogMessage << "* MPI ranks      : "<<GridCmdVectorIntToString(mpi)<<std::endl;
    std::cout<<GridLogMessage << "* Intranode      : "<<GridCmdVectorIntToString(internal)<<std::endl;
    std::cout<<GridLogMessage << "* nodes          : "<<GridCmdVectorIntToString(nodes)<<std::endl;
    std::cout<<GridLogMessage << "* Using "<<threads<<" threads"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;


    ///////// Lattice Init ////////////
    GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

    
    ///////// RNG Init ////////////
    std::vector<int> seeds4({1,2,3,4});
    std::vector<int> seeds5({5,6,7,8});
    GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
    GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
    std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

    ///////// Source preparation ////////////
    LatticeFermion src   (FGrid); random(RNG5,src);
    LatticeFermion ref   (FGrid);
    LatticeFermion tmp   (FGrid);

    RealD N2 = 1.0/::sqrt(norm2(src));
    src = src*N2;
    
    LatticeGaugeField Umu(UGrid);  SU3::HotConfiguration(RNG4,Umu); 

    DomainWallFermionR Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

    ////////////////////////////////////
    // Naive wilson implementation
    ////////////////////////////////////
    {
      LatticeGaugeField Umu5d(FGrid); 
      std::vector<LatticeColourMatrix> U(4,FGrid);
      for(int ss=0;ss<Umu._grid->oSites();ss++){
	for(int s=0;s<Ls;s++){
	  Umu5d._odata[Ls*ss+s] = Umu._odata[ss];
	}
      }
      ref = zero;
      for(int mu=0;mu<Nd;mu++){
	U[mu] = PeekIndex<LorentzIndex>(Umu5d,mu);
      }
      for(int mu=0;mu<Nd;mu++){
	
	tmp = U[mu]*Cshift(src,mu+1,1);
	ref=ref + tmp - Gamma(Gmu[mu])*tmp;
	
	tmp =adj(U[mu])*src;
	tmp =Cshift(tmp,mu+1,-1);
	ref=ref + tmp + Gamma(Gmu[mu])*tmp;
      }
      ref = -0.5*ref;
    }

    LatticeFermion src_e (FrbGrid);
    LatticeFermion src_o (FrbGrid);
    LatticeFermion r_e   (FrbGrid);
    LatticeFermion r_o   (FrbGrid);
    LatticeFermion r_eo  (FGrid);
    LatticeFermion err   (FGrid);
    {

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd,src_o,src);

#if defined(AVX512) 
      const int num_cases = 6;
      std::string fmt("A/S ; A/O ; U/S ; U/O ; G/S ; G/O ");
#else
      const int num_cases = 4;
      std::string fmt("U/S ; U/O ; G/S ; G/O ");
#endif
      controls Cases [] = {
#ifdef AVX512
	{ QCD::WilsonKernelsStatic::OptInlineAsm , QCD::WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptInlineAsm , QCD::WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  },
#endif
	{ QCD::WilsonKernelsStatic::OptHandUnroll, QCD::WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptHandUnroll, QCD::WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptGeneric   , QCD::WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{ QCD::WilsonKernelsStatic::OptGeneric   , QCD::WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  }
      }; 

      for(int c=0;c<num_cases;c++) {

	QCD::WilsonKernelsStatic::Comms = Cases[c].CommsOverlap;
	QCD::WilsonKernelsStatic::Opt   = Cases[c].Opt;
	CartesianCommunicator::SetCommunicatorPolicy(Cases[c].CommsAsynch);

	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
	if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
	if ( sizeof(Real)==4 )   std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
	if ( sizeof(Real)==8 )   std::cout << GridLogMessage<< "* DOUBLE precision "<<std::endl;
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

	int nwarm = 200;
	double t0=usecond();
	FGrid->Barrier();
	for(int i=0;i<nwarm;i++){
	  Dw.DhopEO(src_o,r_e,DaggerNo);
	}
	FGrid->Barrier();
	double t1=usecond();
	//	uint64_t ncall = (uint64_t) 2.5*1000.0*1000.0*nwarm/(t1-t0);
	//	if (ncall < 500) ncall = 500;
	uint64_t ncall = 1000;

	FGrid->Broadcast(0,&ncall,sizeof(ncall));

	//	std::cout << GridLogMessage << " Estimate " << ncall << " calls per second"<<std::endl;
	Dw.ZeroCounters();

	time_statistics timestat;
	std::vector<double> t_time(ncall);
	for(uint64_t i=0;i<ncall;i++){
	  t0=usecond();
	  Dw.DhopEO(src_o,r_e,DaggerNo);
	  t1=usecond();
	  t_time[i] = t1-t0;
	}
	FGrid->Barrier();
	
	double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
	double flops=(1344.0*volume)/2;
	double mf_hi, mf_lo, mf_err;

	timestat.statistics(t_time);
	mf_hi = flops/timestat.min;
	mf_lo = flops/timestat.max;
	mf_err= flops/timestat.min * timestat.err/timestat.mean;

	mflops = flops/timestat.mean;
	mflops_all.push_back(mflops);
	if ( mflops_best == 0   ) mflops_best = mflops;
	if ( mflops_worst== 0   ) mflops_worst= mflops;
	if ( mflops>mflops_best ) mflops_best = mflops;
	if ( mflops<mflops_worst) mflops_worst= mflops;

	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"Deo mflop/s =   "<< mflops << " ("<<mf_err<<") " << mf_lo<<"-"<<mf_hi <<std::endl;
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"Deo mflop/s per rank   "<< mflops/NP<<std::endl;
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"Deo mflop/s per node   "<< mflops/NN<<std::endl;

	Dw.Report();

	Dw.DhopEO(src_o,r_e,DaggerNo);
	Dw.DhopOE(src_e,r_o,DaggerNo);
	setCheckerboard(r_eo,r_o);
	setCheckerboard(r_eo,r_e);
	err = r_eo-ref; 
	std::cout<<GridLogMessage << "norm diff   "<< norm2(err)<<std::endl;
	assert((norm2(err)<1.0e-4));

      }
      robust = mflops_worst/mflops_best;
      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
      std::cout<<GridLogMessage << L<<"^4 x "<<Ls<< " Deo Best  mflop/s        =   "<< mflops_best << " ; " << mflops_best/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage << L<<"^4 x "<<Ls<< " Deo Worst mflop/s        =   "<< mflops_worst<< " ; " << mflops_worst/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage << std::fixed<<std::setprecision(3)<< L<<"^4 x "<<Ls<< " Performance Robustness   =   "<< robust  <<std::endl;
      std::cout<<GridLogMessage <<fmt << std::endl;
      std::cout<<GridLogMessage ;

      for(int i=0;i<mflops_all.size();i++){
	std::cout<<mflops_all[i]/NN<<" ; " ;
      }
      std::cout<<std::endl;
      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    }
    return mflops_best;
  }

};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  CartesianCommunicator::SetCommunicatorPolicy(CartesianCommunicator::CommunicatorPolicySequential);
#ifdef KNL
  LebesgueOrder::Block = std::vector<int>({8,2,2,2});
#else
  LebesgueOrder::Block = std::vector<int>({2,2,2,2});
#endif
  Benchmark::Decomposition();

  int do_memory=1;
  int do_comms =1;
  int do_su3   =0;
  int do_wilson=1;
  int do_dwf   =1;

  if ( do_su3 ) {
    // empty for now
  }
#if 1
  int sel=2;
  std::vector<int> L_list({8,12,16,24});
#else
  int sel=1;
  std::vector<int> L_list({8,12});
#endif
  int selm1=sel-1;
  std::vector<double> robust_list;

  std::vector<double> wilson;
  std::vector<double> dwf4;
  std::vector<double> dwf5;

  if ( do_wilson ) {
    int Ls=1;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Wilson dslash 4D vectorised" <<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    for(int l=0;l<L_list.size();l++){
      double robust;
      wilson.push_back(Benchmark::DWF(1,L_list[l],robust));
    }
  }

  int Ls=16;
  if ( do_dwf ) {
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Domain wall dslash 4D vectorised" <<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    for(int l=0;l<L_list.size();l++){
      double robust;
      double result = Benchmark::DWF(Ls,L_list[l],robust) ;
      dwf4.push_back(result);
      robust_list.push_back(robust);
    }
  }

  if ( do_dwf ) {
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Domain wall dslash 4D vectorised" <<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    for(int l=0;l<L_list.size();l++){
      dwf5.push_back(Benchmark::DWF5(Ls,L_list[l]));
    }

  }

  if ( do_dwf ) {

  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Summary table Ls="<<Ls <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "L \t\t Wilson \t DWF4 \t DWF5 " <<std::endl;
  for(int l=0;l<L_list.size();l++){
    std::cout<<GridLogMessage << L_list[l] <<" \t\t "<< wilson[l]<<" \t "<<dwf4[l]<<" \t "<<dwf5[l] <<std::endl;
  }
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  }

  int NN=NN_global;
  if ( do_memory ) {
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Memory benchmark " <<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    Benchmark::Memory();
  }

  if ( do_comms && (NN>1) ) {
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Communications benchmark " <<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    Benchmark::Comms();
  }

  if ( do_dwf ) {
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Per Node Summary table Ls="<<Ls <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " L \t\t Wilson\t\t DWF4  \t\t DWF5 " <<std::endl;
  for(int l=0;l<L_list.size();l++){
    std::cout<<GridLogMessage << L_list[l] <<" \t\t "<< wilson[l]/NN<<" \t "<<dwf4[l]/NN<<" \t "<<dwf5[l] /NN<<std::endl;
  }
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Comparison point     result: "  << 0.5*(dwf4[sel]+dwf4[selm1])/NN << " Mflop/s per node"<<std::endl;
  std::cout<<GridLogMessage << " Comparison point is 0.5*("<<dwf4[sel]/NN<<"+"<<dwf4[selm1]/NN << ") "<<std::endl;
  std::cout<<std::setprecision(3);
  std::cout<<GridLogMessage << " Comparison point robustness: "  << robust_list[sel] <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

  }


  Grid_finalize();
}
