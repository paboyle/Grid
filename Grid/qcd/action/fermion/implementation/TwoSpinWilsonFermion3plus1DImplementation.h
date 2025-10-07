/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/TwoSpinWilsonFermion2plus1D.cc

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
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/WilsonFermion5D.h>
#include <Grid/perfmon/PerfCount.h>

NAMESPACE_BEGIN(Grid);
  
  // 5d lattice for DWF.
template<class Impl>
TwoSpinWilsonFermion3plus15D<Impl>::TwoSpinWilsonFermion3plus1D(GaugeField &_Umu,
								GridCartesian         &FourDimGrid,
								GridRedBlackCartesian &FourDimRedBlackGrid,
								GridCartesian         &ThreeDimGrid,
								GridRedBlackCartesian &ThreeDimRedBlackGrid,
               RealD _M5,const ImplParams &p) :
  Kernels(p),
  _FourDimGrid        (&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  _ThreeDimGrid        (&ThreeDimGrid),
  _ThreeDimRedBlackGrid(&ThreeDimRedBlackGrid),
  Stencil    (_FourDimGrid,npoint,Even,directions,displacements,p),
  StencilEven(_FourDimRedBlackGrid,npoint,Even,directions,displacements,p), // source is Even
  StencilOdd (_FourDimRedBlackGrid,npoint,Odd ,directions,displacements,p), // source is Odd
  M5(_M5),
  Umu(_ThreeDimGrid),
  UmuEven(_ThreeDimRedBlackGrid),
  UmuOdd (_ThreeDimRedBlackGrid),
  _tmp(&FourDimRedBlackGrid),
  Dirichlet(0)
{
  // some assertions
  assert(FourDimGrid._ndimension==Nd+1);
  assert(ThreeDimGrid._ndimension==Nd);
  assert(ThreeDimRedBlackGrid._ndimension==Nd);
  assert(FourDimRedBlackGrid._ndimension==Nd+1);
  assert(FourDimRedBlackGrid._checker_dim==1); // Don't checker the s direction

  // extent of fifth dim and not spread out
  Ls=FourDimGrid._fdimensions[0];
  assert(FourDimRedBlackGrid._fdimensions[0]==Ls);
  assert(FourDimGrid._processors[0]         ==1);
  assert(FourDimRedBlackGrid._processors[0] ==1);

  // Other dimensions must match the decomposition of the four-D fields 
  for(int d=0;d<Nd;d++){

    assert(FourDimGrid._processors[d+1]         ==ThreeDimGrid._processors[d]);
    assert(FourDimRedBlackGrid._processors[d+1] ==ThreeDimGrid._processors[d]);
    assert(ThreeDimRedBlackGrid._processors[d]   ==ThreeDimGrid._processors[d]);

    assert(FourDimGrid._fdimensions[d+1]        ==ThreeDimGrid._fdimensions[d]);
    assert(FourDimRedBlackGrid._fdimensions[d+1]==ThreeDimGrid._fdimensions[d]);
    assert(ThreeDimRedBlackGrid._fdimensions[d]  ==ThreeDimGrid._fdimensions[d]);

    assert(FourDimGrid._simd_layout[d+1]        ==ThreeDimGrid._simd_layout[d]);
    assert(FourDimRedBlackGrid._simd_layout[d+1]==ThreeDimGrid._simd_layout[d]);
    assert(ThreeDimRedBlackGrid._simd_layout[d]  ==ThreeDimGrid._simd_layout[d]);
  }

  if ( p.dirichlet.size() == Nd+1) {
    Coordinate block = p.dirichlet;
    for(int d=0;d<Nd+1;d++) {
      if ( block[d] ){
	Dirichlet = 1;
	std::cout << GridLogMessage << " WilsonFermion: non-trivial Dirichlet condition "<< block << std::endl;
	std::cout << GridLogMessage << " WilsonFermion: partial Dirichlet "<< p.partialDirichlet << std::endl;
	Block = block;
      }
    }
  } else {
    Coordinate block(Nd+1,0);
    Block = block;
  }

  // Dimension zero of the five-d is the Ls direction
  assert(FourDimRedBlackGrid._simd_layout[0]==1);
  assert(FourDimGrid._simd_layout[0]        ==1);
    
  // Allocate the required comms buffer
  ImportGauge(_Umu);
  // Build lists of exterior only nodes
  int LLs = FourDimGrid._rdimensions[0];
  int vol3;
  vol3=ThreeDimGrid.oSites();
  Stencil.BuildSurfaceList(LLs,vol3);

  vol3=ThreeDimRedBlackGrid.oSites();
  StencilEven.BuildSurfaceList(LLs,vol3);
   StencilOdd.BuildSurfaceList(LLs,vol3);

}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::ImportGauge(const GaugeField &_Umu)
{
  GaugeField HUmu(_Umu.Grid());
  HUmu = _Umu*(-0.5);
  Impl::DoubleStore(GaugeGrid(),Umu,HUmu);
  pickCheckerboard(Even,UmuEven,Umu);
  pickCheckerboard(Odd ,UmuOdd,Umu);
}
template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopDir(const FermionField &in, FermionField &out,int dir5,int disp)
{
  int dir = dir5-1; // Maps to the ordering above in "directions" that is passed to stencil
                    // we drop off the innermost fifth dimension
  //  assert( (disp==1)||(disp==-1) );
  //  assert( (dir>=0)&&(dir<4) ); //must do x,y,z or t;

  int skip = (disp==1) ? 0 : 1;
  int dirdisp = dir+skip*Nd;
  int gamma   = dir+(1-skip)*Nd;

  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in,compressor);
  
  uint64_t Nsite = Umu.Grid()->oSites();
  Kernels::DhopDirKernel(Stencil,Umu,Stencil.CommBuf(),Ls,Nsite,in,out,dirdisp,gamma);

};
template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopDirAll(const FermionField &in, std::vector<FermionField> &out)
{
  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in,compressor);
  uint64_t Nsite = Umu.Grid()->oSites();
  Kernels::DhopDirAll(Stencil,Umu,Stencil.CommBuf(),Ls,Nsite,in,out);
};


template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DerivInternal(StencilImpl & st,
					  DoubledGaugeField & U,
					  GaugeField &mat,
					  const FermionField &A,
					  const FermionField &B,
					  int dag)
{
  assert((dag==DaggerNo) ||(dag==DaggerYes));

  conformable(st.Grid(),A.Grid());
  conformable(st.Grid(),B.Grid());

  Compressor compressor(dag);
  
  FermionField Btilde(B.Grid());
  FermionField Atilde(B.Grid());

  st.HaloExchange(B,compressor);

  Atilde=A;
  int LLs = B.Grid()->_rdimensions[0];


  for (int mu = 0; mu < Nd; mu++) {
    ////////////////////////////////////////////////////////////////////////
    // Flip gamma if dag
    ////////////////////////////////////////////////////////////////////////
    int gamma = mu;
    if (!dag) gamma += Nd;

    ////////////////////////
    // Call the single hop
    ////////////////////////

    int Usites = U.Grid()->oSites();

    Kernels::DhopDirKernel(st, U, st.CommBuf(), Ls, Usites, B, Btilde, mu,gamma);

    ////////////////////////////
    // spin trace outer product
    ////////////////////////////
    Impl::InsertForce5D(mat, Btilde, Atilde, mu);
  }
}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopDeriv(GaugeField &mat,
                                      const FermionField &A,
                                      const FermionField &B,
                                      int dag)
{
  conformable(A.Grid(),FermionGrid());  
  conformable(A.Grid(),B.Grid());

  //conformable(GaugeGrid(),mat.Grid());// this is not general! leaving as a comment

  mat.Checkerboard() = A.Checkerboard();
  //  mat.checkerboard = A.checkerboard;

  DerivInternal(Stencil,Umu,mat,A,B,dag);
}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopDerivEO(GaugeField &mat,
                                        const FermionField &A,
                                        const FermionField &B,
                                        int dag)
{
  conformable(A.Grid(),FermionRedBlackGrid());
  conformable(A.Grid(),B.Grid());

  assert(B.Checkerboard()==Odd);
  assert(A.Checkerboard()==Even);
  mat.Checkerboard() = Even;

  DerivInternal(StencilOdd,UmuEven,mat,A,B,dag);
}


template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopDerivOE(GaugeField &mat,
                                        const FermionField &A,
                                        const FermionField &B,
                                        int dag)
{
  conformable(A.Grid(),FermionRedBlackGrid());
  conformable(A.Grid(),B.Grid());

  assert(B.Checkerboard()==Even);
  assert(A.Checkerboard()==Odd);
  mat.Checkerboard() = Odd;

  DerivInternal(StencilEven,UmuOdd,mat,A,B,dag);
}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopInternal(StencilImpl & st,
                                         DoubledGaugeField & U,
                                         const FermionField &in, FermionField &out,int dag)
{
  DhopInternalSerialComms(st,U,in,out,dag);
}


template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopInternalOverlappedComms(StencilImpl & st,
							DoubledGaugeField & U,
							const FermionField &in, FermionField &out,int dag)
{
  GRID_TRACE("DhopInternalOverlappedComms");
  Compressor compressor(dag);

  int LLs = in.Grid()->_rdimensions[0];
  int len =  U.Grid()->oSites();

  /////////////////////////////
  // Start comms  // Gather intranode and extra node differentiated??
  /////////////////////////////
  {
    //    std::cout << " TwoSpinWilsonFermion3plus1D gather " <<std::endl;
    GRID_TRACE("Gather");
    st.HaloExchangeOptGather(in,compressor); // Put the barrier in the routine
  }
  
  //  std::cout << " TwoSpinWilsonFermion3plus1D Communicate Begin " <<std::endl;
  std::vector<std::vector<CommsRequest_t> > requests;

#if 1
  /////////////////////////////
  // Overlap with comms
  /////////////////////////////
  st.CommunicateBegin(requests);
  st.CommsMergeSHM(compressor);// Could do this inside parallel region overlapped with comms 
#endif

  /////////////////////////////
  // do the compute interior
  /////////////////////////////
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDagInterior");
    Kernels::DhopDagKernel(st,U,st.CommBuf(),LLs,U.oSites(),in,out,1,0);
  } else {
    GRID_TRACE("DhopInterior");
    Kernels::DhopKernel   (st,U,st.CommBuf(),LLs,U.oSites(),in,out,1,0);
  }
  
  //ifdef GRID_ACCELERATED
#if 0
  /////////////////////////////
  // Overlap with comms -- on GPU the interior kernel call is nonblocking
  /////////////////////////////
  st.CommunicateBegin(requests);
  st.CommsMergeSHM(compressor);// Could do this inside parallel region overlapped with comms
#endif
  
  
  /////////////////////////////
  // Complete comms
  /////////////////////////////
  //  std::cout << " TwoSpinWilsonFermion3plus1D Comms Complete " <<std::endl;
  st.CommunicateComplete(requests);
  //  traceStop(id);

  /////////////////////////////
  // do the compute exterior
  /////////////////////////////
  {
    //    std::cout << " TwoSpinWilsonFermion3plus1D Comms Merge " <<std::endl;
    GRID_TRACE("Merge");
    st.CommsMerge(compressor);
  }
  

  //  std::cout << " TwoSpinWilsonFermion3plus1D Exterior " <<std::endl;
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDagExterior");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out,0,1);
  } else {
    GRID_TRACE("DhopExterior");
    Kernels::DhopKernel   (Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out,0,1);
  }
  //  std::cout << " TwoSpinWilsonFermion3plus1D Done " <<std::endl;
}


template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopInternalSerialComms(StencilImpl & st, 
						    DoubledGaugeField & U,
						    const FermionField &in, 
						    FermionField &out,int dag)
{
  GRID_TRACE("DhopInternalSerialComms");
  Compressor compressor(dag);

  int LLs = in.Grid()->_rdimensions[0];

  //  std::cout << " TwoSpinWilsonFermion3plus1D Halo exch " <<std::endl;
  {
    GRID_TRACE("HaloExchange");
    st.HaloExchangeOpt(in,compressor);
  }
  
  //  std::cout << " TwoSpinWilsonFermion3plus1D Dhop " <<std::endl;
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDag");
    Kernels::DhopDagKernel(st,U,st.CommBuf(),LLs,U.oSites(),in,out);
  } else {
    GRID_TRACE("Dhop");
    Kernels::DhopKernel(st,U,st.CommBuf(),LLs,U.oSites(),in,out);
  }
  //  std::cout << " TwoSpinWilsonFermion3plus1D Done " <<std::endl;
}


template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  assert(in.Checkerboard()==Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven,UmuOdd,in,out,dag);
}
template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  assert(in.Checkerboard()==Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd,UmuEven,in,out,dag);
}
template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopComms(const FermionField &in, FermionField &out)
{
  int dag =0 ;
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());
  out.Checkerboard() = in.Checkerboard();
  Compressor compressor(dag);
  Stencil.HaloExchangeOpt(in,compressor);
}
template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DhopCalc(const FermionField &in, FermionField &out,uint64_t *ids)
{
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());

  out.Checkerboard() = in.Checkerboard();

  int LLs = in.Grid()->_rdimensions[0];
  Kernels::DhopKernel(Stencil,Umu,Stencil.CommBuf(),LLs,Umu.oSites(),in,out,ids);
}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil,Umu,in,out,dag);
}
template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::DW(const FermionField &in, FermionField &out,int dag)
{
  out.Checkerboard()=in.Checkerboard();
  Dhop(in,out,dag); // -0.5 is included
  axpy(out,Nd*1.0-M5,in,out);
}
template <class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::Meooe(const FermionField &in, FermionField &out)
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}

template <class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::MeooeDag(const FermionField &in, FermionField &out)
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(Nd*1.0 + M5);
  out = scal * in;
}

template <class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  out = (1.0/(Nd*1.0 + M5))*in;
}

template<class Impl>
void TwoSpinWilsonFermion3plus1D<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in,out);
}
  
NAMESPACE_END(Grid);




