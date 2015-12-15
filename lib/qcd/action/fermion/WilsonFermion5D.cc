#include <Grid.h>
#include <PerfCount.h>

namespace Grid {
namespace QCD {
  
// S-direction is INNERMOST and takes no part in the parity.
const std::vector<int> WilsonFermion5DStatic::directions   ({1,2,3,4, 1, 2, 3, 4});
const std::vector<int> WilsonFermion5DStatic::displacements({1,1,1,1,-1,-1,-1,-1});
int WilsonFermion5DStatic::HandOptDslash;
int WilsonFermion5DStatic::AsmOptDslash;

  // 5d lattice for DWF.
template<class Impl>
WilsonFermion5D<Impl>::WilsonFermion5D(GaugeField &_Umu,
				       GridCartesian         &FiveDimGrid,
				       GridRedBlackCartesian &FiveDimRedBlackGrid,
				       GridCartesian         &FourDimGrid,
				       GridRedBlackCartesian &FourDimRedBlackGrid,
				       RealD _M5,const ImplParams &p) :
  Kernels(p),
  _FiveDimGrid(&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid(&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (_FiveDimGrid,npoint,Even,directions,displacements),
  StencilEven(_FiveDimRedBlackGrid,npoint,Even,directions,displacements), // source is Even
  StencilOdd (_FiveDimRedBlackGrid,npoint,Odd ,directions,displacements), // source is Odd
  M5(_M5),
  Umu(_FourDimGrid),
  UmuEven(_FourDimRedBlackGrid),
  UmuOdd (_FourDimRedBlackGrid),
  Lebesgue(_FourDimGrid),
  LebesgueEvenOdd(_FourDimRedBlackGrid)
{
  // some assertions
  assert(FiveDimGrid._ndimension==5);
  assert(FourDimGrid._ndimension==4);
  
  assert(FiveDimRedBlackGrid._ndimension==5);
  assert(FourDimRedBlackGrid._ndimension==4);

  assert(FiveDimRedBlackGrid._checker_dim==1);

  // Dimension zero of the five-d is the Ls direction
  Ls=FiveDimGrid._fdimensions[0];
  assert(FiveDimRedBlackGrid._fdimensions[0]==Ls);
  assert(FiveDimRedBlackGrid._processors[0] ==1);
  assert(FiveDimRedBlackGrid._simd_layout[0]==1);
  assert(FiveDimGrid._processors[0]         ==1);
  assert(FiveDimGrid._simd_layout[0]        ==1);

  // Other dimensions must match the decomposition of the four-D fields 
  for(int d=0;d<4;d++){
    assert(FourDimRedBlackGrid._fdimensions[d]  ==FourDimGrid._fdimensions[d]);
    assert(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);

    assert(FourDimRedBlackGrid._processors[d]   ==FourDimGrid._processors[d]);
    assert(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);

    assert(FourDimRedBlackGrid._simd_layout[d]  ==FourDimGrid._simd_layout[d]);
    assert(FiveDimRedBlackGrid._simd_layout[d+1]==FourDimGrid._simd_layout[d]);

    assert(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
    assert(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
    assert(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
  }

  // Allocate the required comms buffer
  comm_buf.resize(Stencil._unified_buffer_size); // this is always big enough to contain EO

  ImportGauge(_Umu);
  commtime=0;
  dslashtime=0;
}  
template<class Impl>
void WilsonFermion5D<Impl>::ImportGauge(const GaugeField &_Umu)
{
  Impl::DoubleStore(GaugeGrid(),Umu,_Umu);
  pickCheckerboard(Even,UmuEven,Umu);
  pickCheckerboard(Odd ,UmuOdd,Umu);
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopDir(const FermionField &in, FermionField &out,int dir5,int disp)
{
  int dir = dir5-1; // Maps to the ordering above in "directions" that is passed to stencil
                    // we drop off the innermost fifth dimension
  //  assert( (disp==1)||(disp==-1) );
  //  assert( (dir>=0)&&(dir<4) ); //must do x,y,z or t;

  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in,comm_buf,compressor);
  
  int skip = (disp==1) ? 0 : 1;

  int dirdisp = dir+skip*4;
  int gamma   = dir+(1-skip)*4;

  assert(dirdisp<=7);
  assert(dirdisp>=0);

PARALLEL_FOR_LOOP
  for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      int sU=ss;
      int sF = s+Ls*sU; 
      Kernels::DiracOptDhopDir(Stencil,Umu,comm_buf,sF,sU,in,out,dirdisp,gamma);
    }
  }
};

template<class Impl>
void WilsonFermion5D<Impl>::DerivInternal(StencilImpl & st,
					  DoubledGaugeField & U,
					  GaugeField &mat,
					  const FermionField &A,
					  const FermionField &B,
					  int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));

  conformable(st._grid,A._grid);
  conformable(st._grid,B._grid);

  Compressor compressor(dag);
  
  FermionField Btilde(B._grid);
  FermionField Atilde(B._grid);

  st.HaloExchange(B,comm_buf,compressor);

  Atilde=A;

  for(int mu=0;mu<Nd;mu++){
      
    ////////////////////////////////////////////////////////////////////////
    // Flip gamma if dag
    ////////////////////////////////////////////////////////////////////////
    int gamma = mu;
    if ( !dag ) gamma+= Nd;

    ////////////////////////
    // Call the single hop
    ////////////////////////

PARALLEL_FOR_LOOP
    for(int sss=0;sss<U._grid->oSites();sss++){
      for(int s=0;s<Ls;s++){
	int sU=sss;
	int sF = s+Ls*sU;

	assert ( sF< B._grid->oSites());
	assert ( sU< U._grid->oSites());

	Kernels::DiracOptDhopDir(st,U,comm_buf,sF,sU,B,Btilde,mu,gamma);

    ////////////////////////////
    // spin trace outer product
    ////////////////////////////

      }

    }

    Impl::InsertForce5D(mat,Btilde,Atilde,mu);

  }
}

template<class Impl>
void WilsonFermion5D<Impl>::DhopDeriv(      GaugeField &mat,
					    const FermionField &A,
					    const FermionField &B,
					    int dag)
{
  conformable(A._grid,FermionGrid());  
  conformable(A._grid,B._grid);
  conformable(GaugeGrid(),mat._grid);

  mat.checkerboard = A.checkerboard;

  DerivInternal(Stencil,Umu,mat,A,B,dag);
}

template<class Impl>
void WilsonFermion5D<Impl>::DhopDerivEO(GaugeField &mat,
					const FermionField &A,
					const FermionField &B,
					int dag)
{
  conformable(A._grid,FermionRedBlackGrid());
  conformable(GaugeRedBlackGrid(),mat._grid);
  conformable(A._grid,B._grid);

  assert(B.checkerboard==Odd);
  assert(A.checkerboard==Even);
  mat.checkerboard = Even;

  DerivInternal(StencilOdd,UmuEven,mat,A,B,dag);
}


template<class Impl>
void WilsonFermion5D<Impl>::Report(void)
{
  std::cout<<GridLogMessage << "********************"<<std::endl;
  std::cout<<GridLogMessage << "Halo   time "<<commtime <<" us"<<std::endl;
  std::cout<<GridLogMessage << "Dslash time "<<dslashtime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil All    time "<<Stencil.halotime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "********************"<<std::endl;
  std::cout<<GridLogMessage << "Stencil nosplice time "<<Stencil.nosplicetime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil gather time "<<Stencil.gathertime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil comm   time "<<Stencil.commtime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil scattertime "<<Stencil.scattertime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "********************"<<std::endl;
  std::cout<<GridLogMessage << "Stencil splice time "<<Stencil.splicetime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil comm   time "<<Stencil.commstime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil gathremtime "<<Stencil.gathermtime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil merge  time "<<Stencil.mergetime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "Stencil buf    time "<<Stencil.buftime<<" us"<<std::endl;
  std::cout<<GridLogMessage << "********************"<<std::endl;
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopDerivOE(GaugeField &mat,
				  const FermionField &A,
				  const FermionField &B,
				  int dag)
{
  conformable(A._grid,FermionRedBlackGrid());
  conformable(GaugeRedBlackGrid(),mat._grid);
  conformable(A._grid,B._grid);

  assert(B.checkerboard==Even);
  assert(A.checkerboard==Odd);
  mat.checkerboard = Odd;

  DerivInternal(StencilEven,UmuOdd,mat,A,B,dag);
}

template<class Impl>
void WilsonFermion5D<Impl>::DhopInternal(StencilImpl & st, LebesgueOrder &lo,
					 DoubledGaugeField & U,
					 const FermionField &in, FermionField &out,int dag)
{
  //  assert((dag==DaggerNo) ||(dag==DaggerYes));

  Compressor compressor(dag);

  // Assume balanced KMP_AFFINITY; this is forced in GridThread.h

  int threads = GridThread::GetThreads();
  int HT      = GridThread::GetHyperThreads();
  int cores   = GridThread::GetCores();
  int nwork = U._grid->oSites();
  
  commtime -=usecond();
  st.HaloExchange(in,comm_buf,compressor);
  commtime +=usecond();
  
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  // Not loop ordering and data layout.
  // Designed to create 
  // - per thread reuse in L1 cache for U
  // - 8 linear access unit stride streams per thread for Fermion for hw prefetchable.
  dslashtime -=usecond();
  if ( dag == DaggerYes ) {
    if( this->HandOptDslash ) {
#pragma omp parallel for schedule(static)
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=ss;
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU;
	  Kernels::DiracOptHandDhopSiteDag(st,U,comm_buf,sF,sU,in,out);
	  }
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	{
	  int sd;
	  for(sd=0;sd<Ls;sd++){
	    int sU=ss;
	    int sF = sd+Ls*sU;
	    Kernels::DiracOptDhopSiteDag(st,U,comm_buf,sF,sU,in,out);
	  }
	}
      }
    }
  } else {
    if( this->AsmOptDslash ) {
      //      for(int i=0;i<1;i++){
      //      for(int i=0;i< PerformanceCounter::NumTypes(); i++ ){
      //	PerformanceCounter Counter(i);
      //	Counter.Start();

#pragma omp parallel for 
      for(int t=0;t<threads;t++){

	int hyperthread = t%HT;
	int core        = t/HT;

        int sswork, swork,soff,ssoff,  sU,sF;
	
	GridThread::GetWork(nwork,core,sswork,ssoff,cores);
	GridThread::GetWork(Ls   , hyperthread, swork, soff,HT);

	for(int ss=0;ss<sswork;ss++){
	  for(int s=soff;s<soff+swork;s++){

	    sU=ss+ ssoff;

	    if ( LebesgueOrder::UseLebesgueOrder ) {
	      sU = lo.Reorder(sU);
	    }
	    sF = s+Ls*sU;
	    Kernels::DiracOptAsmDhopSite(st,U,comm_buf,sF,sU,in,out,(uint64_t *)0);// &buf[0]
	  }
	}
      }
      //      Counter.Stop();
      //      Counter.Report();
      //      }
    } else if( this->HandOptDslash ) {
      /*

#pragma omp parallel for schedule(static)
      for(int t=0;t<threads;t++){

	int hyperthread = t%HT;
	int core        = t/HT;

        int sswork, swork,soff,ssoff,  sU,sF;
	
	GridThread::GetWork(nwork,core,sswork,ssoff,cores);
	GridThread::GetWork(Ls   , hyperthread, swork, soff,HT);

	for(int ss=0;ss<sswork;ss++){
	  sU=ss+ ssoff;
	  for(int s=soff;s<soff+swork;s++){
	    sF = s+Ls*sU;
	    Kernels::DiracOptHandDhopSite(st,U,comm_buf,sF,sU,in,out);
	  }
	}
      }
      */

#pragma omp parallel for schedule(static)
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=ss;
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU;
	  Kernels::DiracOptHandDhopSite(st,U,comm_buf,sF,sU,in,out);
	}
      }
    } else { 
PARALLEL_FOR_LOOP
      for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=ss;
	for(int s=0;s<Ls;s++){
	  int sF = s+Ls*sU; 
	  Kernels::DiracOptDhopSite(st,U,comm_buf,sF,sU,in,out);
	}
      }
    }
  }
  dslashtime +=usecond();
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,LebesgueEvenOdd,UmuOdd,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,LebesgueEvenOdd,UmuEven,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  conformable(in._grid,FermionGrid()); // verifies full grid
  conformable(in._grid,out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil,Lebesgue,Umu,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::DW(const FermionField &in, FermionField &out,int dag)
{
  out.checkerboard=in.checkerboard;
  Dhop(in,out,dag); // -0.5 is included
  axpy(out,4.0-M5,in,out);
}

FermOpTemplateInstantiate(WilsonFermion5D);

}}



