/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_trace.h

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
#ifndef GRID_LATTICE_TRACE_H
#define GRID_LATTICE_TRACE_H

///////////////////////////////////////////////
// Tracing, transposing, peeking, poking
///////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace
////////////////////////////////////////////////////////////////////////////////////////////////////
/*
template<class vobj>
inline auto trace(const Lattice<vobj> &lhs)  -> Lattice<decltype(trace(vobj()))>
{
  Lattice<decltype(trace(vobj()))> ret(lhs.Grid());
  autoView(ret_v , ret, AcceleratorWrite);
  autoView(lhs_v , lhs, AcceleratorRead);
  accelerator_for( ss, lhs_v.size(), vobj::Nsimd(), {
    coalescedWrite(ret_v[ss], trace(lhs_v(ss)));
  });
  return ret;
};
*/
    
////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace Index level dependent operation
////////////////////////////////////////////////////////////////////////////////////////////////////
template<int Index,class vobj>
inline auto TraceIndex(const Lattice<vobj> &lhs) -> Lattice<decltype(traceIndex<Index>(vobj()))>
{
  Lattice<decltype(traceIndex<Index>(vobj()))> ret(lhs.Grid());
  autoView( ret_v , ret, AcceleratorWrite);
  autoView( lhs_v , lhs, AcceleratorRead);
  accelerator_for( ss, lhs_v.size(), vobj::Nsimd(), {
    coalescedWrite(ret_v[ss], traceIndex<Index>(lhs_v(ss)));
  });
  return ret;
};

template<int N, class Vec>
Lattice<iScalar<iScalar<iScalar<Vec> > > > Determinant(const Lattice<iScalar<iScalar<iMatrix<Vec, N> > > > &Umu)
{
  GridBase *grid=Umu.Grid();
  auto lvol = grid->lSites();
  Lattice<iScalar<iScalar<iScalar<Vec> > > > ret(grid);
  typedef typename Vec::scalar_type scalar;
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,lvol,{
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);
    Coordinate lcoor;
    grid->LocalIndexToLocalCoor(site, lcoor);
    iScalar<iScalar<iMatrix<scalar, N> > > Us;
    peekLocalSite(Us, Umu_v, lcoor);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	scalar tmp= Us()()(i,j);
	ComplexD ztmp(real(tmp),imag(tmp));
	EigenU(i,j)=ztmp;
      }}
    ComplexD detD  = EigenU.determinant();
    typename Vec::scalar_type det(detD.real(),detD.imag());
    pokeLocalSite(det,ret_v,lcoor);
  });
  return ret;
}

template<int N, class Vec>
Lattice<iScalar<iScalar<iScalar<Vec> > > > Determinant_real(const Lattice<iScalar<iScalar<iMatrix<Vec, N> > > > &Umu)
{
  GridBase *grid=Umu.Grid();
  auto osites = grid->oSites();
  const int Nsimd=grid->Nsimd();
  Lattice<iScalar<iScalar<iScalar<Vec> > > > ret(grid);
  ret.Checkerboard() = Umu.Checkerboard();
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,osites,{
      Eigen::MatrixXd EigenU = Eigen::MatrixXd::Zero(N,N);
      
      iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
      iScalar<iScalar<iScalar<ComplexD> > > Ud;
      for(int lane=0;lane<Nsimd;lane++){
	Us = extractLane(lane,Umu_v[site]);
	for(int i=0;i<N;i++){
	  for(int j=0;j<N;j++){
	    EigenU(i,j)=real(Us()()(i,j));
	  }}
	Ud()()() = EigenU.determinant();
	insertLane(lane,ret_v[site],Ud);
      }
    });
  return ret;
}

template<int N>
Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > Inverse(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
{
#if 0
  GridBase *grid=Umu.Grid();
  auto lvol = grid->lSites();
  Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > ret(grid);
  
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,lvol,{
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);
    Coordinate lcoor;
    grid->LocalIndexToLocalCoor(site, lcoor);
    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;
    peekLocalSite(Us, Umu_v, lcoor);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	EigenU(i,j) = Us()()(i,j);
      }}
    Eigen::MatrixXcd EigenUinv = EigenU.inverse();
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	Ui()()(i,j) = EigenUinv(i,j);
      }}
    pokeLocalSite(Ui,ret_v,lcoor);
  });
#else
  GridBase *grid=Umu.Grid();
  auto osites = grid->oSites();
  const int Nsimd=grid->Nsimd();
  Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > ret(grid);
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,osites,{
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);

    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;

    for(int lane=0;lane<Nsimd;lane++){
      Us = extractLane(lane,Umu_v[site]);
      for(int i=0;i<N;i++){
	for(int j=0;j<N;j++){
	  EigenU(i,j) = Us()()(i,j);
	}}
      Eigen::MatrixXcd EigenUinv = EigenU.inverse();
      for(int i=0;i<N;i++){
	for(int j=0;j<N;j++){
	  Ui()()(i,j) = EigenUinv(i,j);
	}}
      insertLane(lane,ret_v[site],Ui);
    }
  });
#endif  
  return ret;
}

#if 1
/* Helper functions for inversion of real matrix on GPU based on Nobu's code*/
template<class type1,class type2,int N>
accelerator_inline void LUdcmp( iMatrix<type1,N> &LU, iVector<type2,N> &P)
{

  const RealD TINY=1.0e-40;
  
  type1 vv[N], tmp, _max;
  for(int i=0; i<N; i++) vv[i]=0.0;
  
  for(int i=0; i<N; i++){
    _max=0.0;
    for(int j=0; j<N; j++) if( (tmp=abs(LU(i,j))) > _max ) _max = tmp;
    assert( abs(_max) > TINY );
    vv[i] = abs(1.0/_max);
  }

  int imax;
  for(int k=0; k<N; k++){
    _max=0.0;
    for(int i=k; i<N; i++){
      tmp = vv[i] * abs( LU(i,k) );
      if(tmp>_max) {
	_max = tmp;
	imax = i;
      }
    }
    if(k!=imax){
      for(int j=0; j<N; j++){
	tmp = LU(imax,j);
        LU(imax,j) = LU(k,j);
	LU(k,j) = tmp;
      }
      vv[imax] = vv[k];
    }
    P(k)=imax;
    
    for(int i=k+1; i<N; i++){
      LU(i,k) = LU(i,k) / LU(k,k);
      tmp = LU(i,k);
      for(int j=k+1; j<N; j++) LU(i,j) = LU(i,j) - tmp * LU(k,j);
    } // end i
  } // end k
};
  
template<class type1,class type2,int N>
accelerator_inline void solve( iVector<type1,N> &x, const iMatrix<type1,N> LU, const iVector<type2,N> P){
  
  type1 sum = 0.0;

  int ii=0;
  for(int i=0; i<N; i++){
    int ip = P(i);
    sum = x(ip);
    x(ip) = x(i);
    if(ii!=0)for(int j=ii-1;j<i;j++) sum = sum - LU(i,j)*x(j);
    else if (abs(sum)>0.0) ii=i+1;
    x(i) = sum;
  }
  for(int i=N-1; i>=0; i--){
    sum = x(i);
    for(int j=i+1; j<N; j++) sum = sum - LU(i,j)*x(j);
    x(i) = sum/LU(i,i);
  }
};
#endif

template<int N>
Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > Inverse_RealPart(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
{
  GridBase *grid=Umu.Grid();
  auto osites = grid->oSites();
  const int Nsimd=grid->Nsimd();
  Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > ret(grid);
#if 0 // CPU version
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,osites,{
    Eigen::MatrixXd EigenU = Eigen::MatrixXd::Zero(N,N);

    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;

    for(int lane=0;lane<Nsimd;lane++){
      Us = extractLane(lane,Umu_v[site]);
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          EigenU(i,j) = real(Us()()(i,j));
        }}
      Eigen::MatrixXd EigenUinv = EigenU.inverse();
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          Ui()()(i,j) = EigenUinv(i,j);
        }}
      insertLane(lane,ret_v[site],Ui);
    }
  });
#else //GPU version
  autoView(Umu_v,Umu,AcceleratorRead);
  autoView(ret_v,ret,AcceleratorWrite);
#if 1   // For small enough matrices
  accelerator_for(ss,grid->oSites(),Nsimd,{
      iMatrix<RealD, N>  LU;
      iVector<Integer, N> P;
      iVector<RealD, N> e;
      // scalar layout won't coalesce
#ifdef GRID_SIMT
      {
	int blane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
      for(int blane=0;blane<Nsimd;blane++) {
#endif
	
	for(int i=0;i<N;i++){
	  for(int j=0;j<N;j++){
	    LU(i,j) = getlane(toReal(TensorRemove(Umu_v(ss)()()(i,j))),blane);
	  }}
	LUdcmp(LU,P);
	for(int j=0; j<N; j++){
	  for(int i=0; i<N; i++) e(i) = (i==j);
	  solve(e,LU,P);
	  for(int i=0; i<N; i++) putlane(ret_v[ss]()()(i,j),(ComplexD) e(i),blane);
	}
      }
    });
#else   // Eigen supports inversion on GPU's only for matrices of size < 5
  iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;
  accelerator_for(ss,grid->oSites(),Nsimd,{
      Eigen::MatrixXd EigenU = Eigen::MatrixXd::Zero(N,N);
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          EigenU(i,j) = real(Umu_v(ss)()()(i,j));
        }}
      //Linker error occurs w/r/t Eigen when combining the below two lines into a one liner
      Eigen::MatrixXd EigenUinv; 
      EigenUinv= EigenU.inverse();
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          Ui()()(i,j) = EigenUinv(i,j);
        }}
      coalescedWrite(ret_v[ss],Ui);
    });
#endif
#endif
 return ret;
}

NAMESPACE_END(Grid);
#endif

