/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 
    Source file: ./lib/lattice/Lattice_reduction.h
    Copyright (C) 2015
Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#pragma once 
#include <Grid/Grid_Eigen_Dense.h>

#ifdef GRID_WARN_SUBOPTIMAL
#warning "Optimisation alert all these reduction loops are NOT threaded "
#endif     

NAMESPACE_BEGIN(Grid);

template<class vobj>
static void sliceMaddMatrix (Lattice<vobj> &R,Eigen::MatrixXcd &aa,const Lattice<vobj> &X,const Lattice<vobj> &Y,int Orthog,RealD scale=1.0) 
{    
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nblock = X.Grid()->GlobalDimensions()[Orthog];

  GridBase *FullGrid  = X.Grid();
  //  GridBase *SliceGrid = makeSubSliceGrid(FullGrid,Orthog);

  //  Lattice<vobj> Xslice(SliceGrid);
  //  Lattice<vobj> Rslice(SliceGrid);

  assert( FullGrid->_simd_layout[Orthog]==1);

  //FIXME package in a convenient iterator
  //Should loop over a plane orthogonal to direction "Orthog"
  int stride=FullGrid->_slice_stride[Orthog];
  int block =FullGrid->_slice_block [Orthog];
  int nblock=FullGrid->_slice_nblock[Orthog];
  int ostride=FullGrid->_ostride[Orthog];
#pragma omp parallel 
  {
    std::vector<vobj> s_x(Nblock);

#pragma omp for collapse(2)
    for(int n=0;n<nblock;n++){
      for(int b=0;b<block;b++){
	int o  = n*stride + b;

	for(int i=0;i<Nblock;i++){
	  s_x[i] = X[o+i*ostride];
	}

	vobj dot;
	for(int i=0;i<Nblock;i++){
	  dot = Y[o+i*ostride];
	  for(int j=0;j<Nblock;j++){
	    dot = dot + s_x[j]*(scale*aa(j,i));
	  }
	  R[o+i*ostride]=dot;
	}
      }}
  }
};

template<class vobj>
static void sliceMulMatrix (Lattice<vobj> &R,Eigen::MatrixXcd &aa,const Lattice<vobj> &X,int Orthog,RealD scale=1.0) 
{    
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  int Nblock = X.Grid()->GlobalDimensions()[Orthog];

  GridBase *FullGrid  = X.Grid();
  //  GridBase *SliceGrid = makeSubSliceGrid(FullGrid,Orthog);
  //  Lattice<vobj> Xslice(SliceGrid);
  //  Lattice<vobj> Rslice(SliceGrid);

  assert( FullGrid->_simd_layout[Orthog]==1);
  //  int nh =  FullGrid->_ndimension;
  //  int nl = SliceGrid->_ndimension;
  //  int nl=1;

  //FIXME package in a convenient iterator
  //Should loop over a plane orthogonal to direction "Orthog"
  int stride=FullGrid->_slice_stride[Orthog];
  int block =FullGrid->_slice_block [Orthog];
  int nblock=FullGrid->_slice_nblock[Orthog];
  int ostride=FullGrid->_ostride[Orthog];
#pragma omp parallel 
  {
    std::vector<vobj> s_x(Nblock);

#pragma omp for collapse(2)
    for(int n=0;n<nblock;n++){
      for(int b=0;b<block;b++){
	int o  = n*stride + b;

	for(int i=0;i<Nblock;i++){
	  s_x[i] = X[o+i*ostride];
	}

	vobj dot;
	for(int i=0;i<Nblock;i++){
	  dot = s_x[0]*(scale*aa(0,i));
	  for(int j=1;j<Nblock;j++){
	    dot = dot + s_x[j]*(scale*aa(j,i));
	  }
	  R[o+i*ostride]=dot;
	}
      }}
  }

};


template<class vobj>
static void sliceInnerProductMatrix(  Eigen::MatrixXcd &mat, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int Orthog) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  GridBase *FullGrid  = lhs.Grid();
  //  GridBase *SliceGrid = makeSubSliceGrid(FullGrid,Orthog);
  
  int Nblock = FullGrid->GlobalDimensions()[Orthog];
  
  //  Lattice<vobj> Lslice(SliceGrid);
  //  Lattice<vobj> Rslice(SliceGrid);
  
  mat = Eigen::MatrixXcd::Zero(Nblock,Nblock);

  assert( FullGrid->_simd_layout[Orthog]==1);
  //  int nh =  FullGrid->_ndimension;
  //  int nl = SliceGrid->_ndimension;
  //  int nl = nh-1;

  //FIXME package in a convenient iterator
  //Should loop over a plane orthogonal to direction "Orthog"
  int stride=FullGrid->_slice_stride[Orthog];
  int block =FullGrid->_slice_block [Orthog];
  int nblock=FullGrid->_slice_nblock[Orthog];
  int ostride=FullGrid->_ostride[Orthog];

  typedef typename vobj::vector_typeD vector_typeD;

#pragma omp parallel 
  {
    std::vector<vobj> Left(Nblock);
    std::vector<vobj> Right(Nblock);
    Eigen::MatrixXcd  mat_thread = Eigen::MatrixXcd::Zero(Nblock,Nblock);

#pragma omp for collapse(2)
    for(int n=0;n<nblock;n++){
      for(int b=0;b<block;b++){

	int o  = n*stride + b;

	for(int i=0;i<Nblock;i++){
	  Left [i] = lhs[o+i*ostride];
	  Right[i] = rhs[o+i*ostride];
	}

	for(int i=0;i<Nblock;i++){
	  for(int j=0;j<Nblock;j++){
	    auto tmp = innerProduct(Left[i],Right[j]);
	    auto rtmp = TensorRemove(tmp);
	    mat_thread(i,j) += Reduce(rtmp);
	  }}
      }}
#pragma omp critical
    {
      mat += mat_thread;
    }  
  }

  for(int i=0;i<Nblock;i++){
    for(int j=0;j<Nblock;j++){
      ComplexD sum = mat(i,j);
      FullGrid->GlobalSum(sum);
      mat(i,j)=sum;
    }}

  return;
}

NAMESPACE_END(Grid);



