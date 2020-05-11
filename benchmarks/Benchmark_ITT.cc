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

using namespace Grid;

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

    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();

    for(int mu=0;mu<Nd;mu++) if (mpi_layout[mu]>1) nmu++;

    std::vector<double> t_time(Nloop);
    time_statistics timestat;

    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking threaded STENCIL halo exchange in "<<nmu<<" dimensions"<<std::endl;
    std::cout<<GridLogMessage << "===================================================================================================="<<std::endl;
    comms_header();

    for(int lat=16;lat<=maxlat;lat+=8){
      //      for(int Ls=8;Ls<=8;Ls*=2){
      { int Ls=12;

	Coordinate latt_size  ({lat*mpi_layout[0],
	      lat*mpi_layout[1],
	      lat*mpi_layout[2],
	      lat*mpi_layout[3]});
	std::cout << GridLogMessage<< latt_size <<std::endl;
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

	  thread_for(dir,8,{
		     
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
	      tbytes= Grid.StencilSendToRecvFrom((void *)&xbuf[dir][0], xmit_to_rank,
						 (void *)&rbuf[dir][0], recv_from_rank,
						 bytes,dir);
	      thread_critical {
		ncomm++;
		dbytes+=tbytes;
	      }
	    }
          });
	  Grid.Barrier();
	  double stop=usecond();
	  t_time[i] = stop-start; // microseconds
	}

	timestat.statistics(t_time);

	dbytes=dbytes*ppn;
	double xbytes    = dbytes*0.5;
	//	double rbytes    = dbytes*0.5;
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

    Coordinate simd_layout = GridDefaultSimd(Nd,vReal::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();

    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "= Benchmarking a*x + y bandwidth"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "  L  "<<"\t\t"<<"bytes"<<"\t\t\t"<<"GB/s"<<"\t\t"<<"Gflop/s"<<"\t\t seconds"<< "\t\tGB/s / node"<<std::endl;
    std::cout<<GridLogMessage << "----------------------------------------------------------"<<std::endl;
  
    //    uint64_t NP;
    uint64_t NN;


  uint64_t lmax=32;
#define NLOOP (100*lmax*lmax*lmax*lmax/lat/lat/lat/lat)

    GridSerialRNG          sRNG;      sRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    for(int lat=8;lat<=lmax;lat+=8){

      Coordinate latt_size  ({lat*mpi_layout[0],lat*mpi_layout[1],lat*mpi_layout[2],lat*mpi_layout[3]});
      int64_t vol= latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];

      GridCartesian     Grid(latt_size,simd_layout,mpi_layout);

      //      NP= Grid.RankCount();
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
	auto x_v = x.View();
	auto y_v = y.View();
	auto z_v = z.View();
        x_v[0]=z_v[0]; // force serial dependency to prevent optimise away
        y_v[4]=z_v[4];
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


  static double DWF(int Ls,int L)
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
    Coordinate mpi = GridDefaultMpi(); assert(mpi.size()==4);
    Coordinate local({L,L,L,L});

    GridCartesian         * TmpGrid   = SpaceTimeGrid::makeFourDimGrid(Coordinate({72,72,72,72}), 
								       GridDefaultSimd(Nd,vComplex::Nsimd()),
								       GridDefaultMpi());
    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
    NN_global=NN;
    uint64_t SHM=NP/NN;

    Coordinate latt4({local[0]*mpi[0],local[1]*mpi[1],local[2]*mpi[2],local[3]*mpi[3]});

    ///////// Welcome message ////////////
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "Benchmark DWF on "<<L<<"^4 local volume "<<std::endl;
    std::cout<<GridLogMessage << "* Global volume  : "<<GridCmdVectorIntToString(latt4)<<std::endl;
    std::cout<<GridLogMessage << "* Ls             : "<<Ls<<std::endl;
    std::cout<<GridLogMessage << "* ranks          : "<<NP  <<std::endl;
    std::cout<<GridLogMessage << "* nodes          : "<<NN  <<std::endl;
    std::cout<<GridLogMessage << "* ranks/node     : "<<SHM <<std::endl;
    std::cout<<GridLogMessage << "* ranks geom     : "<<GridCmdVectorIntToString(mpi)<<std::endl;
    std::cout<<GridLogMessage << "* Using "<<threads<<" threads"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    ///////// Lattice Init ////////////
    GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
    GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

    
    ///////// RNG Init ////////////
    std::vector<int> seeds4({1,2,3,4});
    std::vector<int> seeds5({5,6,7,8});
    GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
    GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
    std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

    typedef DomainWallFermionF Action;
    typedef typename Action::FermionField Fermion;
    typedef LatticeGaugeFieldF Gauge;
    
    ///////// Source preparation ////////////
    Gauge Umu(UGrid);  SU3::HotConfiguration(RNG4,Umu); 
    Fermion src   (FGrid); random(RNG5,src);
    Fermion src_e (FrbGrid);
    Fermion src_o (FrbGrid);
    Fermion r_e   (FrbGrid);
    Fermion r_o   (FrbGrid);
    Fermion r_eo  (FGrid);
    Action Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

    {

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd,src_o,src);

      const int num_cases = 4;
      std::string fmt("G/S/C ; G/O/C ; G/S/S ; G/O/S ");

      controls Cases [] = {
	{  WilsonKernelsStatic::OptGeneric   ,  WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicyConcurrent  },
	{  WilsonKernelsStatic::OptGeneric   ,  WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicyConcurrent  },
	{  WilsonKernelsStatic::OptGeneric   ,  WilsonKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{  WilsonKernelsStatic::OptGeneric   ,  WilsonKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  }
      }; 

      for(int c=0;c<num_cases;c++) {
	
	WilsonKernelsStatic::Comms = Cases[c].CommsOverlap;
	WilsonKernelsStatic::Opt   = Cases[c].Opt;
	CartesianCommunicator::SetCommunicatorPolicy(Cases[c].CommsAsynch);

	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
	if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
	if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential Comms/Compute" <<std::endl;
	std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

	int nwarm = 10;
	double t0=usecond();
	FGrid->Barrier();
	for(int i=0;i<nwarm;i++){
	  Dw.DhopEO(src_o,r_e,DaggerNo);
	}
	FGrid->Barrier();
	double t1=usecond();
	uint64_t ncall = 50;

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

      }

      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
      std::cout<<GridLogMessage << L<<"^4 x "<<Ls<< " Deo Best  mflop/s        =   "<< mflops_best << " ; " << mflops_best/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage << L<<"^4 x "<<Ls<< " Deo Worst mflop/s        =   "<< mflops_worst<< " ; " << mflops_worst/NN<<" per node " <<std::endl;
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


  static double Staggered(int L)
  {
    double mflops;
    double mflops_best = 0;
    double mflops_worst= 0;
    std::vector<double> mflops_all;

    ///////////////////////////////////////////////////////
    // Set/Get the layout & grid size
    ///////////////////////////////////////////////////////
    int threads = GridThread::GetThreads();
    Coordinate mpi = GridDefaultMpi(); assert(mpi.size()==4);
    Coordinate local({L,L,L,L});
    
    GridCartesian         * TmpGrid   = SpaceTimeGrid::makeFourDimGrid(Coordinate({72,72,72,72}), 
								       GridDefaultSimd(Nd,vComplex::Nsimd()),
								       GridDefaultMpi());
    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
    NN_global=NN;
    uint64_t SHM=NP/NN;

    Coordinate latt4({local[0]*mpi[0],local[1]*mpi[1],local[2]*mpi[2],local[3]*mpi[3]});

    ///////// Welcome message ////////////
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << "Benchmark ImprovedStaggered on "<<L<<"^4 local volume "<<std::endl;
    std::cout<<GridLogMessage << "* Global volume  : "<<GridCmdVectorIntToString(latt4)<<std::endl;
    std::cout<<GridLogMessage << "* ranks          : "<<NP  <<std::endl;
    std::cout<<GridLogMessage << "* nodes          : "<<NN  <<std::endl;
    std::cout<<GridLogMessage << "* ranks/node     : "<<SHM <<std::endl;
    std::cout<<GridLogMessage << "* ranks geom     : "<<GridCmdVectorIntToString(mpi)<<std::endl;
    std::cout<<GridLogMessage << "* Using "<<threads<<" threads"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    ///////// Lattice Init ////////////
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);
    
    ///////// RNG Init ////////////
    std::vector<int> seeds4({1,2,3,4});
    GridParallelRNG          RNG4(FGrid);  RNG4.SeedFixedIntegers(seeds4);
    std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

    RealD mass=0.1;
    RealD c1=9.0/8.0;
    RealD c2=-1.0/24.0;
    RealD u0=1.0;

    typedef ImprovedStaggeredFermionF Action;
    typedef typename Action::FermionField Fermion; 
    typedef LatticeGaugeFieldF Gauge;
    
    Gauge Umu(FGrid);  SU3::HotConfiguration(RNG4,Umu); 

    typename Action::ImplParams params;
    Action Ds(Umu,Umu,*FGrid,*FrbGrid,mass,c1,c2,u0,params);

    ///////// Source preparation ////////////
    Fermion src   (FGrid); random(RNG4,src);
    Fermion src_e (FrbGrid);
    Fermion src_o (FrbGrid);
    Fermion r_e   (FrbGrid);
    Fermion r_o   (FrbGrid);
    Fermion r_eo  (FGrid);
  
    {

      pickCheckerboard(Even,src_e,src);
      pickCheckerboard(Odd,src_o,src);
    
      const int num_cases = 4;
      std::string fmt("G/S/C ; G/O/C ; G/S/S ; G/O/S ");
      
      controls Cases [] = {
	{  StaggeredKernelsStatic::OptGeneric   ,  StaggeredKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicyConcurrent  },
	{  StaggeredKernelsStatic::OptGeneric   ,  StaggeredKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicyConcurrent  },
	{  StaggeredKernelsStatic::OptGeneric   ,  StaggeredKernelsStatic::CommsThenCompute ,CartesianCommunicator::CommunicatorPolicySequential  },
	{  StaggeredKernelsStatic::OptGeneric   ,  StaggeredKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicySequential  }
      }; 

      for(int c=0;c<num_cases;c++) {
	
	StaggeredKernelsStatic::Comms = Cases[c].CommsOverlap;
	StaggeredKernelsStatic::Opt   = Cases[c].Opt;
	CartesianCommunicator::SetCommunicatorPolicy(Cases[c].CommsAsynch);
      
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	if ( StaggeredKernelsStatic::Opt == StaggeredKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc StaggeredKernels" <<std::endl;
	if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
	if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential Comms/Compute" <<std::endl;
	std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	
	int nwarm = 10;
	double t0=usecond();
	FGrid->Barrier();
	for(int i=0;i<nwarm;i++){
	  Ds.DhopEO(src_o,r_e,DaggerNo);
	}
	FGrid->Barrier();
	double t1=usecond();
	uint64_t ncall = 500;

	FGrid->Broadcast(0,&ncall,sizeof(ncall));

	//	std::cout << GridLogMessage << " Estimate " << ncall << " calls per second"<<std::endl;
	Ds.ZeroCounters();

	time_statistics timestat;
	std::vector<double> t_time(ncall);
	for(uint64_t i=0;i<ncall;i++){
	  t0=usecond();
	  Ds.DhopEO(src_o,r_e,DaggerNo);
	  t1=usecond();
	  t_time[i] = t1-t0;
	}
	FGrid->Barrier();
	
	double volume=1;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
	double flops=(1146.0*volume)/2;
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
      
      }

      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
      std::cout<<GridLogMessage << L<<"^4  Deo Best  mflop/s        =   "<< mflops_best << " ; " << mflops_best/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage << L<<"^4  Deo Worst mflop/s        =   "<< mflops_worst<< " ; " << mflops_worst/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage <<fmt << std::endl;
      std::cout<<GridLogMessage ;

      for(int i=0;i<mflops_all.size();i++){
	std::cout<<mflops_all[i]/NN<<" ; " ;
      }
      std::cout<<std::endl;
    }
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
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

  int sel=2;
  std::vector<int> L_list({16,24,32});
  int selm1=sel-1;

  std::vector<double> wilson;
  std::vector<double> dwf4;
  std::vector<double> staggered;

  int Ls=1;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Wilson dslash 4D vectorised" <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  for(int l=0;l<L_list.size();l++){
    wilson.push_back(Benchmark::DWF(Ls,L_list[l]));
  }

  Ls=12;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Domain wall dslash 4D vectorised" <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  for(int l=0;l<L_list.size();l++){
    double result = Benchmark::DWF(Ls,L_list[l]) ;
    dwf4.push_back(result);
  }

  /*
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Improved Staggered dslash 4D vectorised" <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  for(int l=0;l<L_list.size();l++){
    double result = Benchmark::Staggered(L_list[l]) ;
    staggered.push_back(result);
  }
  */

  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << " Summary table Ls="<<Ls <<std::endl;
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
  std::cout<<GridLogMessage << "L \t\t Wilson \t\t DWF4 \t\tt Staggered" <<std::endl;
  for(int l=0;l<L_list.size();l++){
    std::cout<<GridLogMessage << L_list[l] <<" \t\t "<< wilson[l]<<" \t\t "<<dwf4[l] <<std::endl;
  }
  std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

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

    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Per Node Summary table Ls="<<Ls <<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " L \t\t Wilson\t\t DWF4  " <<std::endl;
    for(int l=0;l<L_list.size();l++){
      std::cout<<GridLogMessage << L_list[l] <<" \t\t "<< wilson[l]/NN<<" \t "<<dwf4[l]/NN<<std::endl;
    }
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
    std::cout<<GridLogMessage << " Comparison point     result: "  << 0.5*(dwf4[sel]+dwf4[selm1])/NN << " Mflop/s per node"<<std::endl;
    std::cout<<GridLogMessage << " Comparison point is 0.5*("<<dwf4[sel]/NN<<"+"<<dwf4[selm1]/NN << ") "<<std::endl;
    std::cout<<std::setprecision(3);
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

  Grid_finalize();
}
