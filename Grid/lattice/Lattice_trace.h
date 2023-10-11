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

template<int N>
Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > Inverse(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
{
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
  return ret;
}


NAMESPACE_END(Grid);
#endif

