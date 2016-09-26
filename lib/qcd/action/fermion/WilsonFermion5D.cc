/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonFermion5D.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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
#include <Grid.h>
#include <PerfCount.h>

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
  LebesgueEvenOdd(_FourDimRedBlackGrid)
{
  if (Impl::LsVectorised) { 

    int nsimd = Simd::Nsimd();
    
    // some assertions
    assert(FiveDimGrid._ndimension==5);
    assert(FiveDimRedBlackGrid._ndimension==5);
    assert(FiveDimRedBlackGrid._checker_dim==1); // Don't checker the s direction
    assert(FourDimGrid._ndimension==4);

    // Dimension zero of the five-d is the Ls direction
    Ls=FiveDimGrid._fdimensions[0];
    assert(FiveDimGrid._processors[0]         ==1);
    assert(FiveDimGrid._simd_layout[0]        ==nsimd);

    assert(FiveDimRedBlackGrid._fdimensions[0]==Ls);
    assert(FiveDimRedBlackGrid._processors[0] ==1);
    assert(FiveDimRedBlackGrid._simd_layout[0]==nsimd);

    // Other dimensions must match the decomposition of the four-D fields 
    for(int d=0;d<4;d++){
      assert(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);
      assert(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);
      
      assert(FourDimGrid._simd_layout[d]=1);
      assert(FourDimRedBlackGrid._simd_layout[d]=1);
      assert(FiveDimRedBlackGrid._simd_layout[d+1]==1);

      assert(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
      assert(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
      assert(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
    }

  } else {

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
  }
    
  // Allocate the required comms buffer
  ImportGauge(_Umu);
}
  /*
template<class Impl>
WilsonFermion5D<Impl>::WilsonFermion5D(int simd,GaugeField &_Umu,
				       GridCartesian         &FiveDimGrid,
				       GridRedBlackCartesian &FiveDimRedBlackGrid,
				       GridCartesian         &FourDimGrid,
				       RealD _M5,const ImplParams &p) :
{
  int nsimd = Simd::Nsimd();

  // some assertions
  assert(FiveDimGrid._ndimension==5);
  assert(FiveDimRedBlackGrid._ndimension==5);
  assert(FiveDimRedBlackGrid._checker_dim==0); // Checkerboard the s-direction
  assert(FourDimGrid._ndimension==4);

  // Dimension zero of the five-d is the Ls direction
  Ls=FiveDimGrid._fdimensions[0];
  assert(FiveDimGrid._processors[0]         ==1);
  assert(FiveDimGrid._simd_layout[0]        ==nsimd);

  assert(FiveDimRedBlackGrid._fdimensions[0]==Ls);
  assert(FiveDimRedBlackGrid._processors[0] ==1);
  assert(FiveDimRedBlackGrid._simd_layout[0]==nsimd);

  // Other dimensions must match the decomposition of the four-D fields 
  for(int d=0;d<4;d++){
    assert(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);
    assert(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);

    assert(FourDimGrid._simd_layout[d]=1);
    assert(FiveDimRedBlackGrid._simd_layout[d+1]==1);

    assert(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
    assert(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
    assert(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
  }

  {
  }
}  
  */
     
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

PARALLEL_FOR_LOOP
  for(int ss=0;ss<Umu._grid->oSites();ss++){
    for(int s=0;s<Ls;s++){
      int sU=ss;
      int sF = s+Ls*sU; 
      Kernels::DiracOptDhopDir(Stencil,Umu,Stencil.comm_buf,sF,sU,in,out,dirdisp,gamma);
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

  st.HaloExchange(B,compressor);

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

	Kernels::DiracOptDhopDir(st,U,st.comm_buf,sF,sU,B,Btilde,mu,gamma);

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

  int LLs = in._grid->_rdimensions[0];
  
  st.HaloExchange(in,compressor);
  
  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  if ( dag == DaggerYes ) {
PARALLEL_FOR_LOOP
    for(int ss=0;ss<U._grid->oSites();ss++){
	int sU=ss;
	int sF=LLs*sU;
	Kernels::DiracOptDhopSiteDag(st,lo,U,st.comm_buf,sF,sU,LLs,1,in,out);
    }
  } else {
PARALLEL_FOR_LOOP
    for(int ss=0;ss<U._grid->oSites();ss++){
      int sU=ss;
      int sF=LLs*sU;
      Kernels::DiracOptDhopSite(st,lo,U,st.comm_buf,sF,sU,LLs,1,in,out);
    }
  }
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

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHt(FermionField &out,const FermionField &in) 
{
    // what type LatticeComplex 
    GridBase *_grid = _FourDimGrid;
    conformable(_grid,out._grid);

    typedef typename FermionField::vector_type vector_type;
    typedef typename FermionField::scalar_type ScalComplex;
    typedef iSinglet<ScalComplex> Tcomplex;
    typedef Lattice<iSinglet<vector_type> > LatComplex;

    Gamma::GammaMatrix Gmu [] = {
      Gamma::GammaX,
      Gamma::GammaY,
      Gamma::GammaZ,
      Gamma::GammaT
    };

    std::vector<int> latt_size   = _grid->_fdimensions;


    LatComplex    sk(_grid);  sk = zero;
    LatComplex    sk2(_grid); sk2= zero;

    LatComplex    W(_grid); W= zero;
    LatComplex    a(_grid); a= zero;
    LatComplex    cosha(_grid); a= zero;

    LatComplex     one  (_grid);     one = ScalComplex(1.0,0.0);

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

    W = one - M5 - sk2;
    cosha =  (one + W*W + sk) / (W*2.0);
    
    for(int idx=0;idx<_grid->lSites();idx++){
      Tcomplex cc;
      RealD sgn;
      std::vector<int> lcoor(Nd);
      _grid->LocalIndexToLocalCoor(idx,lcoor);
      peekLocalSite(cc,cosha,lcoor);
      sgn= ::fabs(real(cc));
      cc = ScalComplex(::acosh(sgn),0.0);
      pokeLocalSite(cc,a,lcoor);
    }

  std::cout << "Ht mom space num" << norm2(num)<<std::endl;// num ok
  std::cout << "Ht mom space W  " << norm2(W)<<std::endl;// w ok
  std::cout << "Ht mom space a  " << norm2(a)<<std::endl;// a ok
  denom= ( exp(a) ); 
  std::cout << "Ht mom space den" << norm2(denom)<<std::endl;
  denom= ( exp(a) * W  ); 
  std::cout << "Ht mom space den W" << norm2(denom)<<std::endl;
  denom= ( exp(a) * W - one ); 
  std::cout << "Ht mom space den W- 1" << norm2(denom)<<std::endl;
  denom=where(sk==ScalComplex(0,0),one,denom);
  std::cout << "Ht mom space den W- 1" << denom<<std::endl;
  std::cout << "Ht mom space one " << norm2(one)<<std::endl;
  denom= one/denom;
  std::cout << "Ht mom space div " << norm2(denom)<<std::endl;
  out = num*denom;
  std::cout << "Ht mom space out" << norm2(out)<<std::endl;
}

template<class Impl>
void WilsonFermion5D<Impl>::MomentumSpacePropagatorHw(FermionField &out,const FermionField &in) 
{
    // what type LatticeComplex 
    GridBase *_grid = _FourDimGrid;
    conformable(_grid,out._grid);

    typedef typename FermionField::vector_type vector_type;
    typedef typename FermionField::scalar_type ScalComplex;

    typedef Lattice<iSinglet<vector_type> > LatComplex;

    Gamma::GammaMatrix Gmu [] = {
      Gamma::GammaX,
      Gamma::GammaY,
      Gamma::GammaZ,
      Gamma::GammaT
    };

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

    b_k = sk2 - M5;
     
    w_k = sqrt(sk + b_k*b_k);

    std::cout << "Hw M5 " << M5<<std::endl;
    std::cout << "Hw mom space NUM" << norm2(num)<<std::endl;
    std::cout << "Hw mom space BK" << norm2(b_k)<<std::endl;
    std::cout << "Hw mom space BK" << b_k<<std::endl;
    std::cout << "Hw mom space WK" << norm2(w_k)<<std::endl;
    std::cout << "Hw mom space WK" << w_k<<std::endl;

    denom= ( w_k + b_k ) * (2.0 * M5);
    denom=where(sk==ScalComplex(0,0),one,denom);

    std::cout << "Hw mom space DEMON" << norm2(denom)<<std::endl;

    std::cout << "Hw mom space DEMON" << denom<<std::endl;
    
    std::cout << "Hw mom space one" << norm2(one)<<std::endl;

    denom= one/denom;

    std::cout << "Hw mom space div " << norm2(denom)<<std::endl;

    std::cout << "denom " << denom<<std::endl;

    out = num*denom;
    std::cout << "Hw mom space OUT" << norm2(out)<<std::endl;
}


FermOpTemplateInstantiate(WilsonFermion5D);
GparityFermOpTemplateInstantiate(WilsonFermion5D);
  
}}



