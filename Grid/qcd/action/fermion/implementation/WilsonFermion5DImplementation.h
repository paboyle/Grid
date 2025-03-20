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

NAMESPACE_BEGIN(Grid);
  
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
  Stencil    (_FiveDimGrid,npoint,Even,directions,displacements,p),
  StencilEven(_FiveDimRedBlackGrid,npoint,Even,directions,displacements,p), // source is Even
  StencilOdd (_FiveDimRedBlackGrid,npoint,Odd ,directions,displacements,p), // source is Odd
  M5(_M5),
  Umu(_FourDimGrid),
  UmuEven(_FourDimRedBlackGrid),
  UmuOdd (_FourDimRedBlackGrid),
  _tmp(&FiveDimRedBlackGrid),
  Dirichlet(0)
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

  if ( p.dirichlet.size() == Nd+1) {
    Coordinate block = p.dirichlet;
    if ( block[0] || block[1] || block[2] || block[3] || block[4] ){
      Dirichlet = 1;
      std::cout << GridLogMessage << " WilsonFermion: non-trivial Dirichlet condition "<< block << std::endl;
      std::cout << GridLogMessage << " WilsonFermion: partial Dirichlet "<< p.partialDirichlet << std::endl;
      Block = block;
    }
  } else {
    Coordinate block(Nd+1,0);
    Block = block;
  }

  if (Impl::LsVectorised) { 

    int nsimd = Simd::Nsimd();
    
    // Dimension zero of the five-d is the Ls direction
    assert(FiveDimGrid._simd_layout[0]        ==nsimd);
    assert(FiveDimRedBlackGrid._simd_layout[0]==nsimd);

    for(int d=0;d<4;d++){
      assert(FourDimGrid._simd_layout[d]==1);
      assert(FourDimRedBlackGrid._simd_layout[d]==1);
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

}

template<class Impl>
void WilsonFermion5D<Impl>::ImportGauge(const GaugeField &_Umu)
{
  GaugeField HUmu(_Umu.Grid());
  HUmu = _Umu*(-0.5);
  if ( Dirichlet ) {

    if ( this->Params.partialDirichlet ) {
      std::cout << GridLogMessage << " partialDirichlet BCs " <<Block<<std::endl;
    } else {
      std::cout << GridLogMessage << " FULL Dirichlet BCs " <<Block<<std::endl;
    }
    
    std:: cout << GridLogMessage << "Checking block size multiple of rank boundaries for Dirichlet"<<std::endl;
    for(int d=0;d<Nd;d++) {
      int GaugeBlock = Block[d+1];
      int ldim=GaugeGrid()->LocalDimensions()[d];
      if (GaugeBlock) assert( (GaugeBlock%ldim)==0);
    }

    if (!this->Params.partialDirichlet) {
      std::cout << GridLogMessage << " Dirichlet filtering gauge field BCs block " <<Block<<std::endl;
      Coordinate GaugeBlock(Nd);
      for(int d=0;d<Nd;d++) GaugeBlock[d] = Block[d+1];
      DirichletFilter<GaugeField> Filter(GaugeBlock);
      Filter.applyFilter(HUmu);
    } else {
      std::cout << GridLogMessage << " Dirichlet "<< Dirichlet << " NOT filtered gauge field" <<std::endl;
    }
  }
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

  int skip = (disp==1) ? 0 : 1;
  int dirdisp = dir+skip*4;
  int gamma   = dir+(1-skip)*4;

  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in,compressor);
  
  uint64_t Nsite = Umu.Grid()->oSites();
  Kernels::DhopDirKernel(Stencil,Umu,Stencil.CommBuf(),Ls,Nsite,in,out,dirdisp,gamma);

};
template<class Impl>
void WilsonFermion5D<Impl>::DhopDirAll(const FermionField &in, std::vector<FermionField> &out)
{
  Compressor compressor(DaggerNo);
  Stencil.HaloExchange(in,compressor);
  uint64_t Nsite = Umu.Grid()->oSites();
  Kernels::DhopDirAll(Stencil,Umu,Stencil.CommBuf(),Ls,Nsite,in,out);
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
void WilsonFermion5D<Impl>::DhopDeriv(GaugeField &mat,
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
void WilsonFermion5D<Impl>::DhopDerivEO(GaugeField &mat,
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
void WilsonFermion5D<Impl>::DhopDerivOE(GaugeField &mat,
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
void WilsonFermion5D<Impl>::DhopInternal(StencilImpl & st,
                                         DoubledGaugeField & U,
                                         const FermionField &in, FermionField &out,int dag)
{
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,U,in,out,dag);
  else 
    DhopInternalSerialComms(st,U,in,out,dag);
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopInternalOverlappedComms(StencilImpl & st,
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
    //    std::cout << " WilsonFermion5D gather " <<std::endl;
    GRID_TRACE("Gather");
    st.HaloExchangeOptGather(in,compressor); // Put the barrier in the routine
  }
  
  //  std::cout << " WilsonFermion5D Communicate Begin " <<std::endl;
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
  int Opt = WilsonKernelsStatic::Opt; // Why pass this. Kernels should know
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDagInterior");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out,1,0);
  } else {
    GRID_TRACE("DhopInterior");
    Kernels::DhopKernel   (Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out,1,0);
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
  //  std::cout << " WilsonFermion5D Comms Complete " <<std::endl;
  st.CommunicateComplete(requests);
  //  traceStop(id);

  /////////////////////////////
  // do the compute exterior
  /////////////////////////////
  {
    //    std::cout << " WilsonFermion5D Comms Merge " <<std::endl;
    GRID_TRACE("Merge");
    st.CommsMerge(compressor);
  }
  

  //  std::cout << " WilsonFermion5D Exterior " <<std::endl;
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDagExterior");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out,0,1);
  } else {
    GRID_TRACE("DhopExterior");
    Kernels::DhopKernel   (Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out,0,1);
  }
  //  std::cout << " WilsonFermion5D Done " <<std::endl;
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopInternalSerialComms(StencilImpl & st, 
						    DoubledGaugeField & U,
						    const FermionField &in, 
						    FermionField &out,int dag)
{
  GRID_TRACE("DhopInternalSerialComms");
  Compressor compressor(dag);

  int LLs = in.Grid()->_rdimensions[0];

  //  std::cout << " WilsonFermion5D Halo exch " <<std::endl;
  {
    GRID_TRACE("HaloExchange");
    st.HaloExchangeOpt(in,compressor);
  }
  
  //  std::cout << " WilsonFermion5D Dhop " <<std::endl;
  int Opt = WilsonKernelsStatic::Opt;
  if (dag == DaggerYes) {
    GRID_TRACE("DhopDag");
    Kernels::DhopDagKernel(Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out);
  } else {
    GRID_TRACE("Dhop");
    Kernels::DhopKernel(Opt,st,U,st.CommBuf(),LLs,U.oSites(),in,out);
  }
  //  std::cout << " WilsonFermion5D Done " <<std::endl;
}


template<class Impl>
void WilsonFermion5D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  assert(in.Checkerboard()==Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven,UmuOdd,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  assert(in.Checkerboard()==Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd,UmuEven,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopComms(const FermionField &in, FermionField &out)
{
  int dag =0 ;
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());
  out.Checkerboard() = in.Checkerboard();
  Compressor compressor(dag);
  Stencil.HaloExchangeOpt(in,compressor);
}
template<class Impl>
void WilsonFermion5D<Impl>::DhopCalc(const FermionField &in, FermionField &out,uint64_t *ids)
{
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());

  out.Checkerboard() = in.Checkerboard();

  int LLs = in.Grid()->_rdimensions[0];
  int Opt = WilsonKernelsStatic::Opt;
  Kernels::DhopKernel(Opt,Stencil,Umu,Stencil.CommBuf(),LLs,Umu.oSites(),in,out,ids);
}

template<class Impl>
void WilsonFermion5D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil,Umu,in,out,dag);
}
template<class Impl>
void WilsonFermion5D<Impl>::DW(const FermionField &in, FermionField &out,int dag)
{
  out.Checkerboard()=in.Checkerboard();
  Dhop(in,out,dag); // -0.5 is included
  axpy(out,4.0-M5,in,out);
}

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHt_5d(FermionField &out,const FermionField &in, RealD mass,std::vector<double> twist)
{
  // what type LatticeComplex 
  GridBase *_grid = _FourDimGrid;
  GridBase *_5dgrid = _FiveDimGrid;

  conformable(_5dgrid,out.Grid());

  FermionField   PRsource(_5dgrid);
  FermionField   PLsource(_5dgrid);
  FermionField   buf1_4d(_grid);
  FermionField   buf2_4d(_grid);
  FermionField   GR(_5dgrid);
  FermionField   GL(_5dgrid);
  FermionField   bufL_4d(_grid);
  FermionField   bufR_4d(_grid);

  unsigned int Ls = in.Grid()->_rdimensions[0];
  
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

  Coordinate latt_size   = _grid->_fdimensions;

  LatComplex    sk(_grid);  sk = Zero();
  LatComplex    sk2(_grid); sk2= Zero();
  LatComplex    W(_grid); W= Zero();
  LatComplex    one  (_grid); one = ScalComplex(1.0,0.0);
  LatComplex 	cosha(_grid);
  LatComplex 	kmu(_grid);
  LatComplex 	Wea(_grid);
  LatComplex 	Wema(_grid);
  LatComplex 	ea(_grid);
  LatComplex 	ema(_grid);
  LatComplex 	eaLs(_grid);
  LatComplex 	emaLs(_grid);
  LatComplex 	ea2Ls(_grid);
  LatComplex 	ema2Ls(_grid);
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

  ea = (cosha + sqrt(cosha*cosha-one));
  ema= (cosha - sqrt(cosha*cosha-one));
  eaLs = pow(ea,Ls);
  emaLs= pow(ema,Ls);
  ea2Ls = pow(ea,2.0*Ls);
  ema2Ls= pow(ema,2.0*Ls);
  Wea= abs(W) * ea;
  Wema= abs(W) * ema;
  //  a=log(ea);
  
  sinha = 0.5*(ea - ema);
  sinhaLs = 0.5*(eaLs-emaLs);
  coshaLs = 0.5*(eaLs+emaLs);

  A = one / (abs(W) * sinha * 2.0) * one / (sinhaLs * 2.0);
  F = eaLs * (one - Wea + (Wema - one) * mass*mass);
  F = F + emaLs * (Wema - one + (one - Wea) * mass*mass);
  F = F - abs(W) * sinha * 4.0 * mass;

  Bpp =  (A/F) * (ema2Ls - one) * (one - Wema) * (one - mass*mass * one);
  Bmm =  (A/F) * (one - ea2Ls)  * (one - Wea) * (one - mass*mass * one);
  App =  (A/F) * (ema2Ls - one) * ema * (ema - abs(W)) * (one - mass*mass * one);
  Amm =  (A/F) * (one - ea2Ls)  * ea  * (ea  - abs(W)) * (one - mass*mass * one);
  ABpm = (A/F) * abs(W) * sinha * 2.0  * (one + mass * coshaLs * 2.0 + mass*mass * one);

  //P+ source, P- source
  PRsource = (in + g5 * in) * 0.5;
  PLsource = (in - g5 * in) * 0.5;

  //calculate GR, GL
  for(unsigned int ss=1;ss<=Ls;ss++)
  {
    bufR_4d = Zero();
    bufL_4d = Zero();
    for(unsigned int tt=1;tt<=Ls;tt++)
    {
      //possible sign if W<0
      if((ss+tt)%2==1) signW = abs(W)/W;
      else signW = one;

      unsigned int f = (ss > tt) ? ss-tt : tt-ss; //f = abs(ss-tt)
      //GR
      buf1_4d = Zero();
      ExtractSlice(buf1_4d, PRsource, (tt-1), 0);
      //G(s,t)
      bufR_4d = bufR_4d + A * eaLs * pow(ema,f) * signW * buf1_4d + A * emaLs * pow(ea,f) * signW * buf1_4d;
      //A++*exp(a(s+t))
      bufR_4d = bufR_4d + App * pow(ea,ss) * pow(ea,tt) * signW * buf1_4d ;
      //A+-*exp(a(s-t))
      bufR_4d = bufR_4d + ABpm * pow(ea,ss) * pow(ema,tt) * signW * buf1_4d ;
      //A-+*exp(a(-s+t))
      bufR_4d = bufR_4d + ABpm * pow(ema,ss) * pow(ea,tt) * signW * buf1_4d ;
      //A--*exp(a(-s-t))
      bufR_4d = bufR_4d + Amm * pow(ema,ss) * pow(ema,tt) * signW * buf1_4d ;

      //GL
      buf2_4d = Zero();
      ExtractSlice(buf2_4d, PLsource, (tt-1), 0);
      //G(s,t)
      bufL_4d = bufL_4d + A * eaLs * pow(ema,f) * signW * buf2_4d + A * emaLs * pow(ea,f) * signW * buf2_4d;
      //B++*exp(a(s+t))
      bufL_4d = bufL_4d + Bpp * pow(ea,ss) * pow(ea,tt) * signW * buf2_4d ;
      //B+-*exp(a(s-t))
      bufL_4d = bufL_4d + ABpm * pow(ea,ss) * pow(ema,tt) * signW * buf2_4d ;
      //B-+*exp(a(-s+t))
      bufL_4d = bufL_4d + ABpm * pow(ema,ss) * pow(ea,tt) * signW * buf2_4d ;
      //B--*exp(a(-s-t))
      bufL_4d = bufL_4d + Bmm * pow(ema,ss) * pow(ema,tt) * signW * buf2_4d ;
    }
    InsertSlice(bufR_4d, GR, (ss-1), 0);
    InsertSlice(bufL_4d, GL, (ss-1), 0);
  }

//calculate propagator
  for(unsigned int ss=1;ss<=Ls;ss++)
  {
    bufR_4d = Zero();
    bufL_4d = Zero();

    //(i*gamma_mu*sin(p_mu) - W)*(GL*P- source)
    buf1_4d = Zero();
    ExtractSlice(buf1_4d, GL, (ss-1), 0);
    buf2_4d = Zero();
    for(int mu=0;mu<Nd;mu++) {
      LatticeCoordinate(kmu,mu);
      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
      kmu = TwoPiL * kmu + TwoPiL * one * twist[mu];//twisted boundary
      buf2_4d = buf2_4d + sin(kmu)*ci*(Gamma(Gmu[mu])*buf1_4d);
    }
    bufL_4d = buf2_4d - W * buf1_4d;

    //(i*gamma_mu*sin(p_mu) - W)*(GR*P+ source)
    buf1_4d = Zero();
    ExtractSlice(buf1_4d, GR, (ss-1), 0);
    buf2_4d = Zero();
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
  conformable(_grid,out.Grid());
  
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

  Coordinate latt_size   = _grid->_fdimensions;

  
  FermionField   num  (_grid); num  = Zero();

  LatComplex    sk(_grid);  sk = Zero();
  LatComplex    sk2(_grid); sk2= Zero();
  LatComplex    W(_grid); W= Zero();
  LatComplex    a(_grid); a= Zero();
  LatComplex    one  (_grid); one = ScalComplex(1.0,0.0);
  LatComplex denom(_grid); denom= Zero();
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
  // Cosh alpha -> exp(+/- alpha)
  ////////////////////////////////////////////
  cosha =  (one + W*W + sk) / (abs(W)*2.0);

  Wea = abs(W)*(cosha + sqrt(cosha*cosha-one));
  Wema= abs(W)*(cosha - sqrt(cosha*cosha-one));
  
  num   = num + ( one - Wema ) * mass * in;
  denom= ( Wea - one ) + mass*mass * (one - Wema); 
  out = num/denom;
}

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHw(FermionField &out,const FermionField &in,RealD mass,std::vector<double> twist)
{
  std::vector<double> empty_q(Nd,0.0);
  MomentumSpacePropagatorHwQ(out,in,mass,twist,empty_q);
}
template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHwQ(FermionField &out,const FermionField &in,
						       RealD mass,
						       std::vector<double> twist,
						       std::vector<double> qmu)
{
    Gamma::Algebra Gmu [] = {
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::GammaT
    };

    GridBase *_grid = _FourDimGrid;
    conformable(_grid,out.Grid());

    typedef typename FermionField::vector_type vector_type;
    typedef typename FermionField::scalar_type ScalComplex;

    typedef Lattice<iSinglet<vector_type> > LatComplex;
    typedef iSpinMatrix<ScalComplex> SpinMat;


    Coordinate latt_size   = _grid->_fdimensions;

    LatComplex    sk(_grid);  sk = Zero();
    LatComplex    sk2(_grid); sk2= Zero();

    LatComplex    w_k(_grid); w_k= Zero();
    LatComplex    b_k(_grid); b_k= Zero();

    LatComplex     one  (_grid); one = ScalComplex(1.0,0.0);

    FermionField   num  (_grid); num  = Zero();
    LatComplex denom(_grid); denom= Zero();
    LatComplex kmu(_grid); 
    ScalComplex ci(0.0,1.0);

    std::cout<< "Feynman Rule" << "qmu ("<<qmu[0]<<","<<qmu[1]<<","<<qmu[2]<<","<<qmu[3]<<")"<<std::endl;
    
    for(int mu=0;mu<Nd;mu++) {
      
      LatticeCoordinate(kmu,mu);

      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];

      kmu = TwoPiL * kmu;
      kmu = kmu + TwoPiL * one * twist[mu];//momentum for twisted boundary conditions

      sk2 = sk2 + 2.0*sin(kmu*0.5)*sin(kmu*0.5);

      sk = sk + (sin(kmu)+qmu[mu])*(sin(kmu)+qmu[mu]); 

      // Terms for boosted Fermion
      // 1/2 [ -i gamma.(sin p + q )     ]
      //     [ --------------------- + 1 ]
      //     [         wq + b            ]
      //
      // wq = sqrt( (sinp+q)^2 + b^2 )
      //
      
      num = num - (sin(kmu)+qmu[mu])*ci*(Gamma(Gmu[mu])*in);

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

// Helper macro to reverse Simd vector. Fixme: slow, generic implementation.
#define REVERSE_LS(qSite, qSiteRev, Nsimd) \
{ \
    ExtractBuffer<typename SitePropagator::scalar_object> qSiteVec(Nsimd);	\
    extract(qSite, qSiteVec); \
    for (int i = 0; i < Nsimd / 2; ++i) \
    { \
        typename SitePropagator::scalar_object tmp = qSiteVec[i]; \
        qSiteVec[i] = qSiteVec[Nsimd - i - 1]; \
        qSiteVec[Nsimd - i - 1] = tmp; \
    } \
    merge(qSiteRev, qSiteVec); \
}

 ******************************************************************************/



  
NAMESPACE_END(Grid);




