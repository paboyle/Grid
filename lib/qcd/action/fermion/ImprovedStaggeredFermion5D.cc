/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/ImprovedStaggeredFermion5D.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/ImprovedStaggeredFermion5D.h>
#include <Grid/perfmon/PerfCount.h>

namespace Grid {
namespace QCD {
  
// S-direction is INNERMOST and takes no part in the parity.
const std::vector<int> 
ImprovedStaggeredFermion5DStatic::directions({1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4});
const std::vector<int> 
ImprovedStaggeredFermion5DStatic::displacements({1, 1, 1, 1, -1, -1, -1, -1, 3, 3, 3, 3, -3, -3, -3, -3});

  // 5d lattice for DWF.
template<class Impl>
ImprovedStaggeredFermion5D<Impl>::ImprovedStaggeredFermion5D(GridCartesian         &FiveDimGrid,
							     GridRedBlackCartesian &FiveDimRedBlackGrid,
							     GridCartesian         &FourDimGrid,
							     GridRedBlackCartesian &FourDimRedBlackGrid,
							     RealD _mass,
							     RealD _c1,RealD _c2, RealD _u0,
							     const ImplParams &p) :
  Kernels(p),
  _FiveDimGrid        (&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid        (&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (&FiveDimGrid,npoint,Even,directions,displacements),
  StencilEven(&FiveDimRedBlackGrid,npoint,Even,directions,displacements), // source is Even
  StencilOdd (&FiveDimRedBlackGrid,npoint,Odd ,directions,displacements), // source is Odd
  mass(_mass),
  c1(_c1),
  c2(_c2),
  u0(_u0),
  Umu(&FourDimGrid),
  UmuEven(&FourDimRedBlackGrid),
  UmuOdd (&FourDimRedBlackGrid),
  UUUmu(&FourDimGrid),
  UUUmuEven(&FourDimRedBlackGrid),
  UUUmuOdd(&FourDimRedBlackGrid),
  Lebesgue(&FourDimGrid),
  LebesgueEvenOdd(&FourDimRedBlackGrid),
  _tmp(&FiveDimRedBlackGrid)
{

  // some assertions
  assert(FiveDimGrid._ndimension==5);
  assert(FourDimGrid._ndimension==4);
  assert(FourDimRedBlackGrid._ndimension==4);
  assert(FiveDimRedBlackGrid._ndimension==5);
  assert(FiveDimRedBlackGrid._checker_dim==1); // Don't checker the s direction

  // extent of fifth dim and not spread out
  Ls=FiveDimGrid._fdimensions[0];
  assert(FiveDimRedBlackGrid._fdimensions[0]==Ls);
  assert(FiveDimGrid._processors[0]         ==1);
  assert(FiveDimRedBlackGrid._processors[0] ==1);

  // Other dimensions must match the decomposition of the four-D fields 
  for(int d=0;d<4;d++){
    assert(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
    assert(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);
    assert(FourDimRedBlackGrid._processors[d]   ==FourDimGrid._processors[d]);

    assert(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
    assert(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);
    assert(FourDimRedBlackGrid._fdimensions[d]  ==FourDimGrid._fdimensions[d]);

    assert(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
    assert(FiveDimRedBlackGrid._simd_layout[d+1]==FourDimGrid._simd_layout[d]);
    assert(FourDimRedBlackGrid._simd_layout[d]  ==FourDimGrid._simd_layout[d]);
  }

  if (Impl::LsVectorised) { 

    int nsimd = Simd::Nsimd();
    
    // Dimension zero of the five-d is the Ls direction
    assert(FiveDimGrid._simd_layout[0]        ==nsimd);
    assert(FiveDimRedBlackGrid._simd_layout[0]==nsimd);

    for(int d=0;d<4;d++){
      assert(FourDimGrid._simd_layout[d]=1);
      assert(FourDimRedBlackGrid._simd_layout[d]=1);
      assert(FiveDimRedBlackGrid._simd_layout[d+1]==1);
    }

  } else {
    
    // Dimension zero of the five-d is the Ls direction
    assert(FiveDimRedBlackGrid._simd_layout[0]==1);
    assert(FiveDimGrid._simd_layout[0]        ==1);

  }
  int LLs = FiveDimGrid._rdimensions[0];
  int vol4= FourDimGrid.oSites();
  Stencil.BuildSurfaceList(LLs,vol4);

  vol4=FourDimRedBlackGrid.oSites();
  StencilEven.BuildSurfaceList(LLs,vol4);
  StencilOdd.BuildSurfaceList(LLs,vol4);
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::CopyGaugeCheckerboards(void)
{
  pickCheckerboard(Even, UmuEven,  Umu);
  pickCheckerboard(Odd,  UmuOdd ,  Umu);
  pickCheckerboard(Even, UUUmuEven,UUUmu);
  pickCheckerboard(Odd,  UUUmuOdd, UUUmu);
}
template<class Impl>
ImprovedStaggeredFermion5D<Impl>::ImprovedStaggeredFermion5D(GaugeField &_Uthin,GaugeField &_Ufat,
							     GridCartesian         &FiveDimGrid,
							     GridRedBlackCartesian &FiveDimRedBlackGrid,
							     GridCartesian         &FourDimGrid,
							     GridRedBlackCartesian &FourDimRedBlackGrid,
							     RealD _mass,
							     RealD _c1,RealD _c2, RealD _u0,
							     const ImplParams &p) :
  ImprovedStaggeredFermion5D(FiveDimGrid,FiveDimRedBlackGrid,
			     FourDimGrid,FourDimRedBlackGrid,
			     _mass,_c1,_c2,_u0,p)
{
  ImportGauge(_Uthin,_Ufat);
}

///////////////////////////////////////////////////
// For MILC use; pass three link U's and 1 link U
///////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::ImportGaugeSimple(const GaugeField &_Utriple,const GaugeField &_Ufat) 
{
  /////////////////////////////////////////////////////////////////
  // Trivial import; phases and fattening and such like preapplied
  /////////////////////////////////////////////////////////////////
  for (int mu = 0; mu < Nd; mu++) {

    auto U = PeekIndex<LorentzIndex>(_Utriple, mu);
    Impl::InsertGaugeField(UUUmu,U,mu);

    U = adj( Cshift(U, mu, -3));
    Impl::InsertGaugeField(UUUmu,-U,mu+4);

    U = PeekIndex<LorentzIndex>(_Ufat, mu);
    Impl::InsertGaugeField(Umu,U,mu);

    U = adj( Cshift(U, mu, -1));
    Impl::InsertGaugeField(Umu,-U,mu+4);

  }
  CopyGaugeCheckerboards();
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::ImportGaugeSimple(const DoubledGaugeField &_UUU,const DoubledGaugeField &_U) 
{
  /////////////////////////////////////////////////////////////////
  // Trivial import; phases and fattening and such like preapplied
  /////////////////////////////////////////////////////////////////
  Umu   = _U;
  UUUmu = _UUU;
  CopyGaugeCheckerboards();
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::ImportGauge(const GaugeField &_Uthin,const GaugeField &_Ufat)
{
  ////////////////////////////////////////////////////////
  // Double Store should take two fields for Naik and one hop separately.
  ////////////////////////////////////////////////////////
  Impl::DoubleStore(GaugeGrid(), UUUmu, Umu, _Uthin, _Ufat );

  ////////////////////////////////////////////////////////
  // Apply scale factors to get the right fermion Kinetic term
  // Could pass coeffs into the double store to save work.
  // 0.5 ( U p(x+mu) - Udag(x-mu) p(x-mu) ) 
  ////////////////////////////////////////////////////////
  for (int mu = 0; mu < Nd; mu++) {

    auto U = PeekIndex<LorentzIndex>(Umu, mu);
    PokeIndex<LorentzIndex>(Umu, U*( 0.5*c1/u0), mu );
    
    U = PeekIndex<LorentzIndex>(Umu, mu+4);
    PokeIndex<LorentzIndex>(Umu, U*(-0.5*c1/u0), mu+4);

    U = PeekIndex<LorentzIndex>(UUUmu, mu);
    PokeIndex<LorentzIndex>(UUUmu, U*( 0.5*c2/u0/u0/u0), mu );
    
    U = PeekIndex<LorentzIndex>(UUUmu, mu+4);
    PokeIndex<LorentzIndex>(UUUmu, U*(-0.5*c2/u0/u0/u0), mu+4);
  }

  CopyGaugeCheckerboards();
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDir(const FermionField &in, FermionField &out,int dir5,int disp)
{
  int dir = dir5-1; // Maps to the ordering above in "directions" that is passed to stencil
                    // we drop off the innermost fifth dimension

  Compressor compressor;
  Stencil.HaloExchange(in,compressor);

  parallel_for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      int sU=ss;
      int sF = s+Ls*sU; 
      Kernels::DhopDir(Stencil, Umu, UUUmu, Stencil.CommBuf(), sF, sU, in, out, dir, disp);
    }
  }
};

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DerivInternal(StencilImpl & st,
            DoubledGaugeField & U,
            DoubledGaugeField & UUU,
            GaugeField &mat,
            const FermionField &A,
            const FermionField &B,
            int dag)
{
  // No force terms in multi-rhs solver staggered
  assert(0);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDeriv(GaugeField &mat,
				      const FermionField &A,
				      const FermionField &B,
				      int dag)
{
  assert(0);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDerivEO(GaugeField &mat,
					const FermionField &A,
					const FermionField &B,
					int dag)
{
  assert(0);
}


template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopDerivOE(GaugeField &mat,
					const FermionField &A,
					const FermionField &B,
					int dag)
{
  assert(0);
}

/*CHANGE */
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternal(StencilImpl & st, LebesgueOrder &lo,
						    DoubledGaugeField & U,DoubledGaugeField & UUU,
						    const FermionField &in, FermionField &out,int dag)
{
#ifdef GRID_OMP
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,lo,U,UUU,in,out,dag);
  else
#endif
    DhopInternalSerialComms(st,lo,U,UUU,in,out,dag);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternalOverlappedComms(StencilImpl & st, LebesgueOrder &lo,
								   DoubledGaugeField & U,DoubledGaugeField & UUU,
								   const FermionField &in, FermionField &out,int dag)
{
#ifdef GRID_OMP
  //  assert((dag==DaggerNo) ||(dag==DaggerYes));

  Compressor compressor; 

  int LLs = in._grid->_rdimensions[0];
  int len =  U._grid->oSites();

  DhopFaceTime-=usecond();
  st.Prepare();
  st.HaloGather(in,compressor);
  //  st.HaloExchangeOptGather(in,compressor); // Wilson compressor
  st.CommsMergeSHM(compressor);// Could do this inside parallel region overlapped with comms
  DhopFaceTime+=usecond();

  double ctime=0;
  double ptime=0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Ugly explicit thread mapping introduced for OPA reasons.
  //////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp parallel reduction(max:ctime) reduction(max:ptime)
  {
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int ncomms = CartesianCommunicator::nCommThreads;
    if (ncomms == -1) ncomms = 1;
    assert(nthreads > ncomms);
    if (tid >= ncomms) {
      double start = usecond();
      nthreads -= ncomms;
      int ttid  = tid - ncomms;
      int n     = U._grid->oSites(); // 4d vol
      int chunk = n / nthreads;
      int rem   = n % nthreads;
      int myblock, myn;
      if (ttid < rem) {
        myblock = ttid * chunk + ttid;
        myn = chunk+1;
      } else {
        myblock = ttid*chunk + rem;
        myn = chunk;
      }

      // do the compute
      if (dag == DaggerYes) {
        for (int ss = myblock; ss < myblock+myn; ++ss) {
          int sU = ss;
	  // Interior = 1; Exterior = 0; must implement for staggered
          Kernels::DhopSiteDag(st,lo,U,UUU,st.CommBuf(),LLs,sU,in,out,1,0); //<---------
        }
      } else {
        for (int ss = myblock; ss < myblock+myn; ++ss) {
	  // Interior = 1; Exterior = 0;
          int sU = ss;
          Kernels::DhopSite(st,lo,U,UUU,st.CommBuf(),LLs,sU,in,out,1,0); //<------------
        }
      }
        ptime = usecond() - start;
    } else {
      double start = usecond();
      st.CommunicateThreaded();
      ctime = usecond() - start;
    }
  }
  DhopCommTime += ctime;
  DhopComputeTime+=ptime;

  // First to enter, last to leave timing
  st.CollateThreads();

  DhopFaceTime-=usecond();
  st.CommsMerge(compressor);
  DhopFaceTime+=usecond();

  DhopComputeTime2-=usecond();
  if (dag == DaggerYes) {
    int sz=st.surface_list.size();
    parallel_for (int ss = 0; ss < sz; ss++) {
      int sU = st.surface_list[ss];
      Kernels::DhopSiteDag(st,lo,U,UUU,st.CommBuf(),LLs,sU,in,out,0,1); //<----------
    }
  } else {
    int sz=st.surface_list.size();
    parallel_for (int ss = 0; ss < sz; ss++) {
      int sU = st.surface_list[ss];
      Kernels::DhopSite(st,lo,U,UUU,st.CommBuf(),LLs,sU,in,out,0,1);//<----------
    }
  }
  DhopComputeTime2+=usecond();
#else
  assert(0);
#endif

}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternalSerialComms(StencilImpl & st, LebesgueOrder &lo,
						    DoubledGaugeField & U,DoubledGaugeField & UUU,
						    const FermionField &in, FermionField &out,int dag)
{
  Compressor compressor;
  int LLs = in._grid->_rdimensions[0];



 //double t1=usecond();
  DhopTotalTime -= usecond();
  DhopCommTime -= usecond();
  st.HaloExchange(in,compressor);
  DhopCommTime += usecond();
  
  DhopComputeTime -= usecond();
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  if (dag == DaggerYes) {
    parallel_for (int ss = 0; ss < U._grid->oSites(); ss++) {
      int sU=ss;
      Kernels::DhopSiteDag(st, lo, U, UUU, st.CommBuf(), LLs, sU,in, out);
    }
  } else {
    parallel_for (int ss = 0; ss < U._grid->oSites(); ss++) {
      int sU=ss;
      Kernels::DhopSite(st,lo,U,UUU,st.CommBuf(),LLs,sU,in,out);
    }
  }
  DhopComputeTime += usecond();
  DhopTotalTime   += usecond();
 //double t2=usecond();
 //std::cout << __FILE__ << " " << __func__  << " Total Time " << DhopTotalTime << std::endl;
 //std::cout << __FILE__ << " " << __func__  << " Total Time Org " << t2-t1 << std::endl;
 //std::cout << __FILE__ << " " << __func__  << " Comml Time " << DhopCommTime << std::endl;
 //std::cout << __FILE__ << " " << __func__  << " Compute Time " << DhopComputeTime << std::endl;

}
/*CHANGE END*/

/* ORG
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopInternal(StencilImpl & st, LebesgueOrder &lo,
						    DoubledGaugeField & U,DoubledGaugeField & UUU,
						    const FermionField &in, FermionField &out,int dag)
{
  Compressor compressor;
  int LLs = in._grid->_rdimensions[0];



  DhopTotalTime -= usecond();
  DhopCommTime -= usecond();
  st.HaloExchange(in,compressor);
  DhopCommTime += usecond();
  
  DhopComputeTime -= usecond();
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  if (dag == DaggerYes) {
    parallel_for (int ss = 0; ss < U._grid->oSites(); ss++) {
      int sU=ss;
      Kernels::DhopSiteDag(st, lo, U, UUU, st.CommBuf(), LLs, sU,in, out);
    }
  } else {
    parallel_for (int ss = 0; ss < U._grid->oSites(); ss++) {
      int sU=ss;
	Kernels::DhopSite(st,lo,U,UUU,st.CommBuf(),LLs,sU,in,out);
    }
  }
  DhopComputeTime += usecond();
  DhopTotalTime   += usecond();
}
*/


template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  DhopCalls+=1;
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,LebesgueEvenOdd,UmuOdd,UUUmuOdd,in,out,dag);
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  DhopCalls+=1;
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,LebesgueEvenOdd,UmuEven,UUUmuEven,in,out,dag);
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  DhopCalls+=2;
  conformable(in._grid,FermionGrid()); // verifies full grid
  conformable(in._grid,out._grid);

  out.checkerboard = in.checkerboard;

  DhopInternal(Stencil,Lebesgue,Umu,UUUmu,in,out,dag);
}

template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::Report(void) 
{
  std::vector<int> latt = GridDefaultLatt();          
  RealD volume = Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];
  RealD NP = _FourDimGrid->_Nprocessors;
  RealD NN = _FourDimGrid->NodeCount();

  std::cout << GridLogMessage << "#### Dhop calls report " << std::endl;

  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D Number of DhopEO Calls   : " 
	    << DhopCalls   << std::endl;
  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D TotalTime   /Calls       : " 
	    << DhopTotalTime   / DhopCalls << " us" << std::endl;
  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D CommTime    /Calls       : " 
	    << DhopCommTime    / DhopCalls << " us" << std::endl;
  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D ComputeTime/Calls        : " 
	    << DhopComputeTime / DhopCalls << " us" << std::endl;

  // Average the compute time
  _FourDimGrid->GlobalSum(DhopComputeTime);
  DhopComputeTime/=NP;

  RealD mflops = 1154*volume*DhopCalls/DhopComputeTime/2; // 2 for red black counting
  std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per node       : " << mflops/NN << std::endl;
  
  RealD Fullmflops = 1154*volume*DhopCalls/(DhopTotalTime)/2; // 2 for red black counting
  std::cout << GridLogMessage << "Average mflops/s per call (full)         : " << Fullmflops << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per rank (full): " << Fullmflops/NP << std::endl;
  std::cout << GridLogMessage << "Average mflops/s per call per node (full): " << Fullmflops/NN << std::endl;

  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D Stencil"    <<std::endl;  Stencil.Report();
  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D StencilEven"<<std::endl;  StencilEven.Report();
  std::cout << GridLogMessage << "ImprovedStaggeredFermion5D StencilOdd" <<std::endl;  StencilOdd.Report();
}
template<class Impl>
void ImprovedStaggeredFermion5D<Impl>::ZeroCounters(void) 
{
  DhopCalls       = 0;
  DhopTotalTime    = 0;
  DhopCommTime    = 0;
  DhopComputeTime = 0;
  DhopFaceTime    = 0;


  Stencil.ZeroCounters();
  StencilEven.ZeroCounters();
  StencilOdd.ZeroCounters();
}

/////////////////////////////////////////////////////////////////////////
// Implement the general interface. Here we use SAME mass on all slices
/////////////////////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp) {
  DhopDir(in, out, dir, disp);
}
template <class Impl>
RealD ImprovedStaggeredFermion5D<Impl>::M(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerNo);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
RealD ImprovedStaggeredFermion5D<Impl>::Mdag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Dhop(in, out, DaggerYes);
  return axpy_norm(out, mass, in, out);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Meooe(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MeooeDag(const FermionField &in, FermionField &out) {
  if (in.checkerboard == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::Mooee(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MooeeDag(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  Mooee(in, out);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MooeeInv(const FermionField &in, FermionField &out) {
  out.checkerboard = in.checkerboard;
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::MooeeInvDag(const FermionField &in,
                                      FermionField &out) {
  out.checkerboard = in.checkerboard;
  MooeeInv(in, out);
}

//////////////////////////////////////////////////////// 
// Conserved current - not yet implemented.
////////////////////////////////////////////////////////
template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
                                                         PropagatorField &q_in_2,
                                                         PropagatorField &q_out,
                                                         Current curr_type,
                                                         unsigned int mu)
{
    assert(0);
}

template <class Impl>
void ImprovedStaggeredFermion5D<Impl>::SeqConservedCurrent(PropagatorField &q_in, 
                                              PropagatorField &q_out,
                                              Current curr_type,
                                              unsigned int mu,
                                              unsigned int tmin, 
                                              unsigned int tmax,
					      ComplexField &lattice_cmplx)
{
    assert(0);

}

FermOpStaggeredTemplateInstantiate(ImprovedStaggeredFermion5D);
FermOpStaggeredVec5dTemplateInstantiate(ImprovedStaggeredFermion5D);
  
}}



