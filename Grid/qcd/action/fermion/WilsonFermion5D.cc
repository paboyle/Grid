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
Author: Andrew Lawson <andrew.lawson1991@gmail.com>
Author: Vera Guelpers <V.M.Guelpers@soton.ac.uk>

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

   //  std::cout << GridLogMessage << " SurfaceLists "<< Stencil.surface_list.size()
   //                       <<" " << StencilEven.surface_list.size()<<std::endl;

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
  if ( DhopCalls > 0){
    std::cout << GridLogMessage << "WilsonFermion5D Stencil     Reporti()"    <<std::endl;  Stencil.Reporti(DhopCalls);
    std::cout << GridLogMessage << "WilsonFermion5D StencilEven Reporti()"<<std::endl;  StencilEven.Reporti(DhopCalls);
    std::cout << GridLogMessage << "WilsonFermion5D StencilOdd  Reporti()" <<std::endl;  StencilOdd.Reporti(DhopCalls);
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
  Stencil.ZeroCountersi();
  StencilEven.ZeroCountersi();
  StencilOdd.ZeroCountersi();
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
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHt_5d(FermionField &out,const FermionField &in, RealD mass,std::vector<double> twist)
{
  // what type LatticeComplex 
  GridBase *_grid = _FourDimGrid;
  GridBase *_5dgrid = _FiveDimGrid;

  conformable(_5dgrid,out._grid);

  FermionField   PRsource(_5dgrid);
  FermionField   PLsource(_5dgrid);
  FermionField   buf1_4d(_grid);
  FermionField   buf2_4d(_grid);
  FermionField   GR(_5dgrid);
  FermionField   GL(_5dgrid);
  FermionField   bufL_4d(_grid);
  FermionField   bufR_4d(_grid);

  unsigned int Ls = in._grid->_rdimensions[0];
  
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

  Gamma g5(Gamma::Algebra::Gamma5);

  std::vector<int> latt_size   = _grid->_fdimensions;

  LatComplex    sk(_grid);  sk = zero;
  LatComplex    sk2(_grid); sk2= zero;
  LatComplex    W(_grid); W= zero;
  LatComplex    a(_grid); a= zero;
  LatComplex    one  (_grid); one = ScalComplex(1.0,0.0);
  LatComplex 	cosha(_grid);
  LatComplex 	kmu(_grid);
  LatComplex 	Wea(_grid);
  LatComplex 	Wema(_grid);
  LatComplex 	sinha(_grid);
  LatComplex 	sinhaLs(_grid);
  LatComplex 	coshaLs(_grid);
  LatComplex 	A(_grid);
  LatComplex 	F(_grid);
  LatComplex 	App(_grid);
  LatComplex 	Amm(_grid);
  LatComplex 	Bpp(_grid);
  LatComplex 	Bmm(_grid);
  LatComplex 	ABpm(_grid); //Apm=Amp=Bpm=Bmp
  LatComplex 	signW(_grid);

  ScalComplex ci(0.0,1.0);

  for(int mu=0;mu<Nd;mu++) {
    
    LatticeCoordinate(kmu,mu);
    
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    
    kmu = TwoPiL * kmu;
    kmu = kmu + TwoPiL * one * twist[mu];//momentum for twisted boundary conditions
    
    sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
    sk  = sk  +     sin(kmu)    *sin(kmu);
  }
  
  W = one - M5 + sk2;

  ////////////////////////////////////////////
  // Cosh alpha -> alpha
  ////////////////////////////////////////////
  cosha = (one + W*W + sk) / (abs(W)*2.0);

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

  Wea = ( exp( a) * abs(W)  );
  Wema= ( exp(-a) * abs(W)  );
  sinha = 0.5*(exp( a) - exp(-a));
  sinhaLs = 0.5*(exp( a*Ls) - exp(-a*Ls));
  coshaLs = 0.5*(exp( a*Ls) + exp(-a*Ls));

  A = one / (abs(W) * sinha * 2.0) * one / (sinhaLs * 2.0);
  F = exp( a*Ls) * (one - Wea + (Wema - one) * mass*mass);
  F = F + exp(-a*Ls) * (Wema - one + (one - Wea) * mass*mass);
  F = F - abs(W) * sinha * 4.0 * mass;

  Bpp =  (A/F) * (exp(-a*Ls*2.0) - one) * (one - Wema) * (one - mass*mass * one);
  Bmm =  (A/F) * (one - exp(a*Ls*2.0)) * (one - Wea) * (one - mass*mass * one);
  App =  (A/F) * (exp(-a*Ls*2.0) - one) * exp(-a) * (exp(-a) - abs(W)) * (one - mass*mass * one);
  Amm =  (A/F) * (one - exp(a*Ls*2.0)) * exp(a) * (exp(a) - abs(W)) * (one - mass*mass * one);
  ABpm = (A/F) * abs(W) * sinha * 2.0  * (one + mass * coshaLs * 2.0 + mass*mass * one);

  //P+ source, P- source
  PRsource = (in + g5 * in) * 0.5;
  PLsource = (in - g5 * in) * 0.5;

  //calculate GR, GL
  for(unsigned int ss=1;ss<=Ls;ss++)
  {
    bufR_4d = zero;
    bufL_4d = zero;
    for(unsigned int tt=1;tt<=Ls;tt++)
    {
      //possible sign if W<0
      if((ss+tt)%2==1) signW = abs(W)/W;
      else signW = one;

      unsigned int f = (ss > tt) ? ss-tt : tt-ss; //f = abs(ss-tt)
      //GR
      buf1_4d = zero;
      ExtractSlice(buf1_4d, PRsource, (tt-1), 0);
      //G(s,t)
      bufR_4d = bufR_4d + A * exp(a*Ls) * exp(-a*f) * signW * buf1_4d + A * exp(-a*Ls) * exp(a*f) * signW * buf1_4d;
      //A++*exp(a(s+t))
      bufR_4d = bufR_4d + App * exp(a*ss) * exp(a*tt) * signW * buf1_4d ;
      //A+-*exp(a(s-t))
      bufR_4d = bufR_4d + ABpm * exp(a*ss) * exp(-a*tt) * signW * buf1_4d ;
      //A-+*exp(a(-s+t))
      bufR_4d = bufR_4d + ABpm * exp(-a*ss) * exp(a*tt) * signW * buf1_4d ;
      //A--*exp(a(-s-t))
      bufR_4d = bufR_4d + Amm * exp(-a*ss) * exp(-a*tt) * signW * buf1_4d ;

      //GL
      buf2_4d = zero;
      ExtractSlice(buf2_4d, PLsource, (tt-1), 0);
      //G(s,t)
      bufL_4d = bufL_4d + A * exp(a*Ls) * exp(-a*f) * signW * buf2_4d + A * exp(-a*Ls) * exp(a*f) * signW * buf2_4d;
      //B++*exp(a(s+t))
      bufL_4d = bufL_4d + Bpp * exp(a*ss) * exp(a*tt) * signW * buf2_4d ;
      //B+-*exp(a(s-t))
      bufL_4d = bufL_4d + ABpm * exp(a*ss) * exp(-a*tt) * signW * buf2_4d ;
      //B-+*exp(a(-s+t))
      bufL_4d = bufL_4d + ABpm * exp(-a*ss) * exp(a*tt) * signW * buf2_4d ;
      //B--*exp(a(-s-t))
      bufL_4d = bufL_4d + Bmm * exp(-a*ss) * exp(-a*tt) * signW * buf2_4d ;
    }
    InsertSlice(bufR_4d, GR, (ss-1), 0);
    InsertSlice(bufL_4d, GL, (ss-1), 0);
  }

//calculate propagator
  for(unsigned int ss=1;ss<=Ls;ss++)
  {
    bufR_4d = zero;
    bufL_4d = zero;

    //(i*gamma_mu*sin(p_mu) - W)*(GL*P- source)
    buf1_4d = zero;
    ExtractSlice(buf1_4d, GL, (ss-1), 0);
    buf2_4d = zero;
    for(int mu=0;mu<Nd;mu++) {
      LatticeCoordinate(kmu,mu);
      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      kmu = TwoPiL * kmu + TwoPiL * one * twist[mu];//twisted boundary
      buf2_4d = buf2_4d + sin(kmu)*ci*(Gamma(Gmu[mu])*buf1_4d);
    }
    bufL_4d = buf2_4d - W * buf1_4d;

    //(i*gamma_mu*sin(p_mu) - W)*(GR*P+ source)
    buf1_4d = zero;
    ExtractSlice(buf1_4d, GR, (ss-1), 0);
    buf2_4d = zero;
    for(int mu=0;mu<Nd;mu++) {
      LatticeCoordinate(kmu,mu);
      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      kmu = TwoPiL * kmu + TwoPiL * one * twist[mu];//twisted boundary
      buf2_4d = buf2_4d + sin(kmu)*ci*(Gamma(Gmu[mu])*buf1_4d);
    }
    bufR_4d = buf2_4d - W * buf1_4d;

    //(delta(s-1,u) - m*delta(s,1)*delta(u,Ls))*GL
    if(ss==1){
      ExtractSlice(buf1_4d, GL, (Ls-1), 0);
      bufL_4d = bufL_4d - mass*buf1_4d;
    }
    else {
      ExtractSlice(buf1_4d, GL, (ss-2), 0);
      bufL_4d = bufL_4d + buf1_4d;
    }

    //(delta(s+1,u) - m*delta(s,Ls)*delta(u,1))*GR
    if(ss==Ls){
      ExtractSlice(buf1_4d, GR, 0, 0);
      bufR_4d = bufR_4d - mass*buf1_4d;
    }
    else {
      ExtractSlice(buf1_4d, GR, ss, 0);
      bufR_4d = bufR_4d + buf1_4d;
    }
    buf1_4d = bufL_4d + bufR_4d;
    InsertSlice(buf1_4d, out, (ss-1), 0);
  }


  out = out * (-1.0);
}

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHt(FermionField &out,const FermionField &in, RealD mass,std::vector<double> twist)
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
    kmu = kmu + TwoPiL * one * twist[mu];//momentum for twisted boundary conditions
    
    sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);
    sk  = sk  +     sin(kmu)    *sin(kmu); 
    
    num = num - sin(kmu)*ci*(Gamma(Gmu[mu])*in);
    
  }
  
  W = one - M5 + sk2;

  ////////////////////////////////////////////
  // Cosh alpha -> alpha
  ////////////////////////////////////////////
  cosha =  (one + W*W + sk) / (abs(W)*2.0);

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
  
  Wea = ( exp( a) * abs(W)  );
  Wema= ( exp(-a) * abs(W)  );
  
  num   = num + ( one - Wema ) * mass * in;
  denom= ( Wea - one ) + mass*mass * (one - Wema); 
  out = num/denom;
}

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHw(FermionField &out,const FermionField &in,RealD mass,std::vector<double> twist)
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
      kmu = kmu + TwoPiL * one * twist[mu];//momentum for twisted boundary conditions

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

/*******************************************************************************
 * Conserved current utilities for Wilson fermions, for contracting propagators
 * to make a conserved current sink or inserting the conserved current 
 * sequentially.
 ******************************************************************************/

// Helper macro to reverse Simd vector. Fixme: slow, generic implementation.
#define REVERSE_LS(qSite, qSiteRev, Nsimd) \
{ \
    std::vector<typename SitePropagator::scalar_object> qSiteVec(Nsimd); \
    extract(qSite, qSiteVec); \
    for (int i = 0; i < Nsimd / 2; ++i) \
    { \
        typename SitePropagator::scalar_object tmp = qSiteVec[i]; \
        qSiteVec[i] = qSiteVec[Nsimd - i - 1]; \
        qSiteVec[Nsimd - i - 1] = tmp; \
    } \
    merge(qSiteRev, qSiteVec); \
}

//          psi = chiralProjectPlus(Result_s[Ls/2-1]);
//          psi+= chiralProjectMinus(Result_s[Ls/2]);
//         PJ5q+=localInnerProduct(psi,psi);

template<class vobj> 
Lattice<vobj> spProj5p(const Lattice<vobj> & in)
{
  GridBase *grid=in._grid;
  Gamma G5(Gamma::Algebra::Gamma5);
  Lattice<vobj> ret(grid);
  parallel_for(int ss=0;ss<grid->oSites();ss++){
    ret._odata[ss] = in._odata[ss] + G5*in._odata[ss];
  }
  return ret;
}
template<class vobj> 
Lattice<vobj> spProj5m(const Lattice<vobj> & in)
{
  Gamma G5(Gamma::Algebra::Gamma5);
  GridBase *grid=in._grid;
  Lattice<vobj> ret(grid);
  parallel_for(int ss=0;ss<grid->oSites();ss++){
    ret._odata[ss] = in._odata[ss] - G5*in._odata[ss];
  }
  return ret;
}

template <class Impl>
void WilsonFermion5D<Impl>::ContractJ5q(FermionField &q_in,ComplexField &J5q)
{
  conformable(GaugeGrid(), J5q._grid);
  conformable(q_in._grid, FermionGrid());

  // 4d field
  int Ls = this->Ls;
  FermionField psi(GaugeGrid());
  FermionField p_plus (GaugeGrid());
  FermionField p_minus(GaugeGrid());
  FermionField p(GaugeGrid());

  ExtractSlice(p_plus , q_in, Ls/2   , 0);
  ExtractSlice(p_minus, q_in, Ls/2-1 , 0);
  p_plus = spProj5p(p_plus );
  p_minus= spProj5m(p_minus);
  p=p_plus+p_minus;
  J5q = localInnerProduct(p,p);
}

template <class Impl>
void WilsonFermion5D<Impl>::ContractJ5q(PropagatorField &q_in,ComplexField &J5q)
{
  conformable(GaugeGrid(), J5q._grid);
  conformable(q_in._grid, FermionGrid());

  // 4d field
  int Ls = this->Ls;
  PropagatorField psi(GaugeGrid());
  PropagatorField p_plus (GaugeGrid());
  PropagatorField p_minus(GaugeGrid());
  PropagatorField p(GaugeGrid());

  ExtractSlice(p_plus , q_in, Ls/2   , 0);
  ExtractSlice(p_minus, q_in, Ls/2-1 , 0);
  p_plus = spProj5p(p_plus );
  p_minus= spProj5m(p_minus);
  p=p_plus+p_minus;
  J5q = localInnerProduct(p,p);
}

template <class Impl>
void WilsonFermion5D<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
                                                     PropagatorField &q_in_2,
                                                     PropagatorField &q_out,
                                                     Current curr_type,
                                                     unsigned int mu)
{
    conformable(q_in_1._grid, FermionGrid());
    conformable(q_in_1._grid, q_in_2._grid);
    conformable(_FourDimGrid, q_out._grid);

    PropagatorField tmp1(FermionGrid()), tmp2(FermionGrid());
    unsigned int LLs = q_in_1._grid->_rdimensions[0];
    q_out = zero;

    // Forward, need q1(x + mu, s), q2(x, Ls - 1 - s). Backward, need q1(x, s), 
    // q2(x + mu, Ls - 1 - s). 5D lattice so shift 4D coordinate mu by one.
    tmp1 = Cshift(q_in_1, mu + 1, 1);
    tmp2 = Cshift(q_in_2, mu + 1, 1);
    parallel_for (unsigned int sU = 0; sU < Umu._grid->oSites(); ++sU)
    {
        unsigned int sF1 = sU * LLs;
        unsigned int sF2 = (sU + 1) * LLs - 1;

        for (unsigned int s = 0; s < LLs; ++s)
        {
            bool axial_sign = ((curr_type == Current::Axial) && \
                               (s < (LLs / 2)));
            SitePropagator qSite2, qmuSite2;

            // If vectorised in 5th dimension, reverse q2 vector to match up
            // sites correctly.
            if (Impl::LsVectorised)
            {
                REVERSE_LS(q_in_2._odata[sF2], qSite2, Ls / LLs);
                REVERSE_LS(tmp2._odata[sF2], qmuSite2, Ls / LLs);
            }
            else
            {
                qSite2   = q_in_2._odata[sF2];
                qmuSite2 = tmp2._odata[sF2];
            }
            Kernels::ContractConservedCurrentSiteFwd(tmp1._odata[sF1], 
                                                     qSite2, 
                                                     q_out._odata[sU],
                                                     Umu, sU, mu, axial_sign);
            Kernels::ContractConservedCurrentSiteBwd(q_in_1._odata[sF1],
                                                     qmuSite2,
                                                     q_out._odata[sU],
                                                     Umu, sU, mu, axial_sign);
            sF1++;
            sF2--;
        }
    }
}


template <class Impl>
void WilsonFermion5D<Impl>::SeqConservedCurrent(PropagatorField &q_in, 
                                                PropagatorField &q_out,
                                                Current curr_type, 
                                                unsigned int mu,
                                                unsigned int tmin, 
                                                unsigned int tmax,
						ComplexField &lattice_cmplx)
{
    conformable(q_in._grid, FermionGrid());
    conformable(q_in._grid, q_out._grid);
    PropagatorField tmp(GaugeGrid()),tmp2(GaugeGrid());
    unsigned int tshift = (mu == Tp) ? 1 : 0;
    unsigned int LLs = q_in._grid->_rdimensions[0];
    unsigned int LLt    = GridDefaultLatt()[Tp];

    q_out = zero;
    LatticeInteger coords(_FourDimGrid);
    LatticeCoordinate(coords, Tp);


    for (unsigned int s = 0; s < LLs; ++s)
    {
        bool axial_sign = ((curr_type == Current::Axial) && (s < (LLs / 2)));
	bool tadpole_sign = (curr_type == Current::Tadpole);
	bool switch_sgn = tadpole_sign || axial_sign;


        //forward direction: Need q(x + mu, s)*A(x)
        ExtractSlice(tmp2, q_in, s, 0);  //q(x,s) 
        tmp = Cshift(tmp2, mu, 1);	 //q(x+mu,s)
        tmp2 = tmp*lattice_cmplx;	 //q(x+mu,s)*A(x)	

    	parallel_for (unsigned int sU = 0; sU < Umu._grid->oSites(); ++sU)
    	{
            // Compute the sequential conserved current insertion only if our simd
            // object contains a timeslice we need.
            vInteger t_mask   = ((coords._odata[sU] >= tmin) &&
                	         (coords._odata[sU] <= tmax));
            Integer timeSlices = Reduce(t_mask);

            if (timeSlices > 0)
            {
		unsigned int sF = sU * LLs + s;
                Kernels::SeqConservedCurrentSiteFwd(tmp2._odata[sU], 
                                              q_out._odata[sF], Umu, sU,
                                              mu, t_mask, switch_sgn);
            }

        }

        //backward direction: Need q(x - mu, s)*A(x-mu)
        ExtractSlice(tmp2, q_in, s, 0);  //q(x,s)
        tmp = lattice_cmplx*tmp2;	 //q(x,s)*A(x)
        tmp2 = Cshift(tmp, mu, -1);	 //q(x-mu,s)*A(x-mu,s)

    	parallel_for (unsigned int sU = 0; sU < Umu._grid->oSites(); ++sU)
    	{
            vInteger  t_mask     = ((coords._odata[sU] >= (tmin + tshift)) && 
                   	  	    (coords._odata[sU] <= (tmax + tshift)));

	    //if tmax = LLt-1 (last timeslice) include timeslice 0 if the time is shifted (mu=3)	
	    unsigned int t0 = 0;
	    if((tmax==LLt-1) && (tshift==1)) t_mask = (t_mask || (coords._odata[sU] == t0 ));

            Integer timeSlices = Reduce(t_mask);

            if (timeSlices > 0)
            {
		unsigned int sF = sU * LLs + s; 
        	Kernels::SeqConservedCurrentSiteBwd(tmp2._odata[sU], 
                                             q_out._odata[sF], Umu, sU,
                                             mu, t_mask, axial_sign);
            }
	}
    }
}





FermOpTemplateInstantiate(WilsonFermion5D);
GparityFermOpTemplateInstantiate(WilsonFermion5D);
  
}}



