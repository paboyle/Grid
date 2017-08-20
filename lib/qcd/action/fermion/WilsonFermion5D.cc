/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonFermion5D.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
#include <Grid/qcd/action/fermion/WilsonFermion5D.h>
#include <Grid/perfmon/PerfCount.h>

namespace Grid {
namespace QCD {
  
// S-direction is INNERMOST and takes no part in the parity.
const std::vector<int> WilsonFermion5DStatic::directions   ({1,2,3,4, 1, 2, 3, 4});
const std::vector<int> WilsonFermion5DStatic::displacements({1,1,1,1,-1,-1,-1,-1});

  // 5d lattice for DWF.
template<class Impl>
WilsonFermion5D<Impl>::WilsonFermion5D(GaugeField &_Umu,
               GridCartesian         &FiveDimGrid,
               GridRedBlackCartesian &FiveDimRedBlackGrid,
               GridCartesian         &FourDimGrid,
               GridRedBlackCartesian &FourDimRedBlackGrid,
               RealD _M5,const ImplParams &p) :
  Kernels(p),
  _FiveDimGrid        (&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid        (&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (_FiveDimGrid,npoint,Even,directions,displacements),
  StencilEven(_FiveDimRedBlackGrid,npoint,Even,directions,displacements), // source is Even
  StencilOdd (_FiveDimRedBlackGrid,npoint,Odd ,directions,displacements), // source is Odd
  M5(_M5),
  Umu(_FourDimGrid),
  UmuEven(_FourDimRedBlackGrid),
  UmuOdd (_FourDimRedBlackGrid),
  Lebesgue(_FourDimGrid),
  LebesgueEvenOdd(_FourDimRedBlackGrid),
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
    
  // Allocate the required comms buffer
  ImportGauge(_Umu);
  // Build lists of exterior only nodes
  int LLs = FiveDimGrid._rdimensions[0];
  int vol4;
  vol4=FourDimGrid.oSites();
  Stencil.BuildSurfaceList(LLs,vol4);
  vol4=FourDimRedBlackGrid.oSites();
  StencilEven.BuildSurfaceList(LLs,vol4);
   StencilOdd.BuildSurfaceList(LLs,vol4);

  std::cout << GridLogMessage << " SurfaceLists "<< Stencil.surface_list.size()
                       <<" " << StencilEven.surface_list.size()<<std::endl;

}
     
template<class Impl>
void WilsonFermion5D<Impl>::Report(void)
{
  RealD NP     = _FourDimGrid->_Nprocessors;
  RealD NN     = _FourDimGrid->NodeCount();
  RealD volume = Ls;  
  std::vector<int> latt = _FourDimGrid->GlobalDimensions();
  for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];

  if ( DhopCalls > 0 ) {
    std::cout << GridLogMessage << "#### Dhop calls report " << std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D Number of DhopEO Calls   : " << DhopCalls   << std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D TotalTime   /Calls        : " << DhopTotalTime   / DhopCalls << " us" << std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D CommTime    /Calls        : " << DhopCommTime    / DhopCalls << " us" << std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D FaceTime    /Calls        : " << DhopFaceTime    / DhopCalls << " us" << std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D ComputeTime1/Calls        : " << DhopComputeTime / DhopCalls << " us" << std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D ComputeTime2/Calls        : " << DhopComputeTime2/ DhopCalls << " us" << std::endl;

    // Average the compute time
    _FourDimGrid->GlobalSum(DhopComputeTime);
    DhopComputeTime/=NP;
    RealD mflops = 1344*volume*DhopCalls/DhopComputeTime/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank       : " << mflops/NP << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per node       : " << mflops/NN << std::endl;

    RealD Fullmflops = 1344*volume*DhopCalls/(DhopTotalTime)/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call (full)         : " << Fullmflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per rank (full): " << Fullmflops/NP << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per node (full): " << Fullmflops/NN << std::endl;

   }

  if ( DerivCalls > 0 ) {
    std::cout << GridLogMessage << "#### Deriv calls report "<< std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D Number of Deriv Calls    : " <<DerivCalls <<std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D CommTime/Calls           : " <<DerivCommTime/DerivCalls<<" us" <<std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D ComputeTime/Calls        : " <<DerivComputeTime/DerivCalls<<" us" <<std::endl;
    std::cout << GridLogMessage << "WilsonFermion5D Dhop ComputeTime/Calls   : " <<DerivDhopComputeTime/DerivCalls<<" us" <<std::endl;
    
    RealD mflops = 144*volume*DerivCalls/DerivDhopComputeTime;
    std::cout << GridLogMessage << "Average mflops/s per call                : " << mflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per node       : " << mflops/NP << std::endl;

    RealD Fullmflops = 144*volume*DerivCalls/(DerivDhopComputeTime+DerivCommTime)/2; // 2 for red black counting
    std::cout << GridLogMessage << "Average mflops/s per call (full)         : " << Fullmflops << std::endl;
    std::cout << GridLogMessage << "Average mflops/s per call per node (full): " << Fullmflops/NP << std::endl;  }

  if (DerivCalls > 0 || DhopCalls > 0){
    std::cout << GridLogMessage << "WilsonFermion5D Stencil"    <<std::endl;  Stencil.Report();
    std::cout << GridLogMessage << "WilsonFermion5D StencilEven"<<std::endl;  StencilEven.Report();
    std::cout << GridLogMessage << "WilsonFermion5D StencilOdd" <<std::endl;  StencilOdd.Report();
  }
}

template<class Impl>
void WilsonFermion5D<Impl>::ZeroCounters(void) {
  DhopCalls       = 0;
  DhopCommTime    = 0;
  DhopComputeTime = 0;
  DhopComputeTime2= 0;
  DhopFaceTime    = 0;
  DhopTotalTime   = 0;

  DerivCalls       = 0;
  DerivCommTime    = 0;
  DerivComputeTime = 0;
  DerivDhopComputeTime = 0;

  Stencil.ZeroCounters();
  StencilEven.ZeroCounters();
  StencilOdd.ZeroCounters();
}


template<class Impl>
void WilsonFermion5D<Impl>::ImportGauge(const GaugeField &_Umu)
{
  GaugeField HUmu(_Umu._grid);
  HUmu = _Umu*(-0.5);
  Impl::DoubleStore(GaugeGrid(),Umu,HUmu);
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
  Stencil.HaloExchange(in,compressor);
  
  int skip = (disp==1) ? 0 : 1;

  int dirdisp = dir+skip*4;
  int gamma   = dir+(1-skip)*4;

  assert(dirdisp<=7);
  assert(dirdisp>=0);

  parallel_for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      int sU=ss;
      int sF = s+Ls*sU; 
      Kernels::DhopDir(Stencil,Umu,Stencil.CommBuf(),sF,sU,in,out,dirdisp,gamma);
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
  DerivCalls++;
  assert((dag==DaggerNo) ||(dag==DaggerYes));

  conformable(st._grid,A._grid);
  conformable(st._grid,B._grid);

  Compressor compressor(dag);
  
  FermionField Btilde(B._grid);
  FermionField Atilde(B._grid);

  DerivCommTime-=usecond();
  st.HaloExchange(B,compressor);
  DerivCommTime+=usecond();

  Atilde=A;
  int LLs = B._grid->_rdimensions[0];


  DerivComputeTime-=usecond();
  for (int mu = 0; mu < Nd; mu++) {
    ////////////////////////////////////////////////////////////////////////
    // Flip gamma if dag
    ////////////////////////////////////////////////////////////////////////
    int gamma = mu;
    if (!dag) gamma += Nd;

    ////////////////////////
    // Call the single hop
    ////////////////////////

    DerivDhopComputeTime -= usecond();
    parallel_for (int sss = 0; sss < U._grid->oSites(); sss++) {
      for (int s = 0; s < Ls; s++) {
        int sU = sss;
        int sF = s + Ls * sU;

        assert(sF < B._grid->oSites());
        assert(sU < U._grid->oSites());

        Kernels::DhopDir(st, U, st.CommBuf(), sF, sU, B, Btilde, mu, gamma);

        ////////////////////////////
        // spin trace outer product
        ////////////////////////////
      }
    }
    ////////////////////////////
    // spin trace outer product
    ////////////////////////////
    DerivDhopComputeTime += usecond();
    Impl::InsertForce5D(mat, Btilde, Atilde, mu);
  }
  DerivComputeTime += usecond();
}

template<class Impl>
void WilsonFermion5D<Impl>::DhopDeriv(GaugeField &mat,
                                      const FermionField &A,
                                      const FermionField &B,
                                      int dag)
{
  conformable(A._grid,FermionGrid());  
  conformable(A._grid,B._grid);

  //conformable(GaugeGrid(),mat._grid);// this is not general! leaving as a comment

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
  conformable(A._grid,B._grid);

  assert(B.checkerboard==Odd);
  assert(A.checkerboard==Even);
  mat.checkerboard = Even;

  DerivInternal(StencilOdd,UmuEven,mat,A,B,dag);
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopDerivOE(GaugeField &mat,
                                        const FermionField &A,
                                        const FermionField &B,
                                        int dag)
{
  conformable(A._grid,FermionRedBlackGrid());
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
  DhopTotalTime-=usecond();
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,lo,U,in,out,dag);
  else 
#endif
    DhopInternalSerialComms(st,lo,U,in,out,dag);
  DhopTotalTime+=usecond();
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopInternalOverlappedComms(StencilImpl & st, LebesgueOrder &lo,
							DoubledGaugeField & U,
							const FermionField &in, FermionField &out,int dag)
{
#ifdef GRID_OMP
  //  assert((dag==DaggerNo) ||(dag==DaggerYes));

  Compressor compressor(dag);

  int LLs = in._grid->_rdimensions[0];
  int len =  U._grid->oSites();

  DhopFaceTime-=usecond();
  st.HaloExchangeOptGather(in,compressor);
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
      int ttid = tid - ncomms;
      int n = U._grid->oSites();
      int chunk = n / nthreads;
      int rem = n % nthreads;
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
	  int sF = LLs * sU;
	  Kernels::DhopSiteDag(st,lo,U,st.CommBuf(),sF,sU,LLs,1,in,out,1,0);
	}
      } else {
	for (int ss = myblock; ss < myblock+myn; ++ss) {
	  int sU = ss;
	  int sF = LLs * sU;
	  Kernels::DhopSite(st,lo,U,st.CommBuf(),sF,sU,LLs,1,in,out,1,0);
	}
      }
	ptime = usecond() - start;
    }
    {
      double start = usecond();
      st.CommunicateThreaded();
      ctime = usecond() - start;
    }
  }
  DhopCommTime += ctime;
  DhopComputeTime+=ptime;

  DhopFaceTime-=usecond();
  st.CommsMerge(compressor);
  DhopFaceTime+=usecond();

  DhopComputeTime2-=usecond();
  if (dag == DaggerYes) {
    int sz=st.surface_list.size();
    parallel_for (int ss = 0; ss < sz; ss++) {
      int sU = st.surface_list[ss];
      int sF = LLs * sU;
      Kernels::DhopSiteDag(st,lo,U,st.CommBuf(),sF,sU,LLs,1,in,out,0,1);
    }
  } else {
    int sz=st.surface_list.size();
    parallel_for (int ss = 0; ss < sz; ss++) {
      int sU = st.surface_list[ss];
      int sF = LLs * sU;
      Kernels::DhopSite(st,lo,U,st.CommBuf(),sF,sU,LLs,1,in,out,0,1);
    }
  }
  DhopComputeTime2+=usecond();
#else 
  assert(0);
#endif
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopInternalSerialComms(StencilImpl & st, LebesgueOrder &lo,
					 DoubledGaugeField & U,
					 const FermionField &in, FermionField &out,int dag)
{
  //  assert((dag==DaggerNo) ||(dag==DaggerYes));
  Compressor compressor(dag);

  int LLs = in._grid->_rdimensions[0];
  
  DhopCommTime-=usecond();
  st.HaloExchangeOpt(in,compressor);
  DhopCommTime+=usecond();
  
  DhopComputeTime-=usecond();
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion

  if (dag == DaggerYes) {
    parallel_for (int ss = 0; ss < U._grid->oSites(); ss++) {
      int sU = ss;
      int sF = LLs * sU;
      Kernels::DhopSiteDag(st,lo,U,st.CommBuf(),sF,sU,LLs,1,in,out);
    }
  } else {
    parallel_for (int ss = 0; ss < U._grid->oSites(); ss++) {
      int sU = ss;
      int sF = LLs * sU;
      Kernels::DhopSite(st,lo,U,st.CommBuf(),sF,sU,LLs,1,in,out);
    }
  }
  DhopComputeTime+=usecond();
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  DhopCalls++;
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Even);
  out.checkerboard = Odd;

  DhopInternal(StencilEven,LebesgueEvenOdd,UmuOdd,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  DhopCalls++;
  conformable(in._grid,FermionRedBlackGrid());    // verifies half grid
  conformable(in._grid,out._grid); // drops the cb check

  assert(in.checkerboard==Odd);
  out.checkerboard = Even;

  DhopInternal(StencilOdd,LebesgueEvenOdd,UmuEven,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  DhopCalls+=2;
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

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHt(FermionField &out,const FermionField &in, RealD mass) 
{
  // what type LatticeComplex 
  GridBase *_grid = _FourDimGrid;
  conformable(_grid,out._grid);
  
  typedef typename FermionField::vector_type vector_type;
  typedef typename FermionField::scalar_type ScalComplex;
  typedef iSinglet<ScalComplex> Tcomplex;
  typedef Lattice<iSinglet<vector_type> > LatComplex;
  
  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

  std::vector<int> latt_size   = _grid->_fdimensions;

  
  FermionField   num  (_grid); num  = zero;

  LatComplex    sk(_grid);  sk = zero;
  LatComplex    sk2(_grid); sk2= zero;
  LatComplex    W(_grid); W= zero;
  LatComplex    a(_grid); a= zero;
  LatComplex    one  (_grid); one = ScalComplex(1.0,0.0);
  LatComplex denom(_grid); denom= zero;
  LatComplex cosha(_grid); 
  LatComplex kmu(_grid); 
  LatComplex Wea(_grid); 
  LatComplex Wema(_grid); 

  ScalComplex ci(0.0,1.0);
  
  for(int mu=0;mu<Nd;mu++) {
    
    LatticeCoordinate(kmu,mu);
    
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    
    kmu = TwoPiL * kmu;
    
    sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
    sk  = sk  +     sin(kmu)    *sin(kmu); 
    
    num = num - sin(kmu)*ci*(Gamma(Gmu[mu])*in);
    
  }
  
  W = one - M5 + sk2;

  ////////////////////////////////////////////
  // Cosh alpha -> alpha
  ////////////////////////////////////////////
  cosha =  (one + W*W + sk) / (W*2.0);

  // FIXME Need a Lattice acosh
  for(int idx=0;idx<_grid->lSites();idx++){
    std::vector<int> lcoor(Nd);
    Tcomplex cc;
    RealD sgn;
    _grid->LocalIndexToLocalCoor(idx,lcoor);
    peekLocalSite(cc,cosha,lcoor);
    assert((double)real(cc)>=1.0);
    assert(fabs((double)imag(cc))<=1.0e-15);
    cc = ScalComplex(::acosh(real(cc)),0.0);
    pokeLocalSite(cc,a,lcoor);
  }
  
  Wea = ( exp( a) * W  ); 
  Wema= ( exp(-a) * W  ); 
  
  num   = num + ( one - Wema ) * mass * in;
  denom= ( Wea - one ) + mass*mass * (one - Wema); 
  out = num/denom;
}

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHw(FermionField &out,const FermionField &in,RealD mass) 
{
    Gamma::Algebra Gmu [] = {
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::GammaT
    };

    GridBase *_grid = _FourDimGrid;
    conformable(_grid,out._grid);

    typedef typename FermionField::vector_type vector_type;
    typedef typename FermionField::scalar_type ScalComplex;

    typedef Lattice<iSinglet<vector_type> > LatComplex;


    std::vector<int> latt_size   = _grid->_fdimensions;

    LatComplex    sk(_grid);  sk = zero;
    LatComplex    sk2(_grid); sk2= zero;

    LatComplex    w_k(_grid); w_k= zero;
    LatComplex    b_k(_grid); b_k= zero;

    LatComplex     one  (_grid); one = ScalComplex(1.0,0.0);

    FermionField   num  (_grid); num  = zero;
    LatComplex denom(_grid); denom= zero;
    LatComplex kmu(_grid); 
    ScalComplex ci(0.0,1.0);

    for(int mu=0;mu<Nd;mu++) {

      LatticeCoordinate(kmu,mu);

      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];

      kmu = TwoPiL * kmu;

      sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
      sk  = sk  + sin(kmu)*sin(kmu); 

      num = num - sin(kmu)*ci*(Gamma(Gmu[mu])*in);

    }
    num = num + mass * in ;

    b_k = sk2 - M5;
     
    w_k = sqrt(sk + b_k*b_k);

    denom= ( w_k + b_k + mass*mass) ;

    denom= one/denom;
    out = num*denom;

}

FermOpTemplateInstantiate(WilsonFermion5D);
GparityFermOpTemplateInstantiate(WilsonFermion5D);
  
}}



