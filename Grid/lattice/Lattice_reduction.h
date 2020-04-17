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


#ifdef GRID_NVCC
#include <Grid/lattice/Lattice_reduction_gpu.h>
#endif

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////
// FIXME this should promote to double and accumulate
//////////////////////////////////////////////////////
template<class vobj>
inline typename vobj::scalar_object sum_cpu(const vobj *arg, Integer osites)
{
  typedef typename vobj::scalar_object  sobj;

  const int Nsimd = vobj::Nsimd();
  const int nthread = GridThread::GetThreads();

  Vector<sobj> sumarray(nthread);
  for(int i=0;i<nthread;i++){
    sumarray[i]=Zero();
  }
  
  thread_for(thr,nthread, {
    int nwork, mywork, myoff;
    nwork = osites;
    GridThread::GetWork(nwork,thr,mywork,myoff);
    vobj vvsum=Zero();
    for(int ss=myoff;ss<mywork+myoff; ss++){
      vvsum = vvsum + arg[ss];
    }
    sumarray[thr]=Reduce(vvsum);
  });
  
  sobj ssum=Zero();  // sum across threads
  for(int i=0;i<nthread;i++){
    ssum = ssum+sumarray[i];
  } 
  
  return ssum;
}
template<class vobj>
inline typename vobj::scalar_object sum(const vobj *arg, Integer osites)
{
#ifdef GRID_NVCC
  return sum_gpu(arg,osites);
#else
  return sum_cpu(arg,osites);
#endif  
}
template<class vobj>
inline typename vobj::scalar_object sum(const Lattice<vobj> &arg)
{
  auto arg_v = arg.View();
  Integer osites = arg.Grid()->oSites();
  auto ssum= sum(&arg_v[0],osites);
  arg.Grid()->GlobalSum(ssum);
  return ssum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Deterministic Reduction operations
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vobj> inline RealD norm2(const Lattice<vobj> &arg){
  ComplexD nrm = innerProduct(arg,arg);
  return real(nrm); 
}

// Double inner product
template<class vobj>
inline ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right)
{
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_typeD vector_type;
  ComplexD  nrm;
  
  GridBase *grid = left.Grid();
  
  // Might make all code paths go this way.
  auto left_v = left.View();
  auto right_v=right.View();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();
  
#ifdef GRID_NVCC
  // GPU - SIMT lane compliance...
  typedef decltype(innerProduct(left_v[0],right_v[0])) inner_t;
  Vector<inner_t> inner_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
  

  accelerator_for( ss, sites, nsimd,{
      auto x_l = left_v(ss);
      auto y_l = right_v(ss);
      coalescedWrite(inner_tmp_v[ss],innerProduct(x_l,y_l));
  })

  // This is in single precision and fails some tests
  // Need a sumD that sums in double
  nrm = TensorRemove(sumD_gpu(inner_tmp_v,sites));  
#else
  // CPU 
  typedef decltype(innerProductD(left_v[0],right_v[0])) inner_t;
  Vector<inner_t> inner_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
  
  accelerator_for( ss, sites, nsimd,{
      auto x_l = left_v[ss];
      auto y_l = right_v[ss];
      inner_tmp_v[ss]=innerProductD(x_l,y_l);
  })
  nrm = TensorRemove(sum(inner_tmp_v,sites));
#endif
  grid->GlobalSum(nrm);

  return nrm;
}

/////////////////////////
// Fast axpby_norm
// z = a x + b y
// return norm z
/////////////////////////
template<class sobj,class vobj> strong_inline RealD 
axpy_norm_fast(Lattice<vobj> &z,sobj a,const Lattice<vobj> &x,const Lattice<vobj> &y) 
{
  sobj one(1.0);
  return axpby_norm_fast(z,a,one,x,y);
}

template<class sobj,class vobj> strong_inline RealD 
axpby_norm_fast(Lattice<vobj> &z,sobj a,sobj b,const Lattice<vobj> &x,const Lattice<vobj> &y) 
{
  z.Checkerboard() = x.Checkerboard();
  conformable(z,x);
  conformable(x,y);

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_typeD vector_type;
  RealD  nrm;
  
  GridBase *grid = x.Grid();

  auto x_v=x.View();
  auto y_v=y.View();
  auto z_v=z.View();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();
  
#ifdef GRID_NVCC
  // GPU
  typedef decltype(innerProduct(x_v[0],y_v[0])) inner_t;
  Vector<inner_t> inner_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];

  accelerator_for( ss, sites, nsimd,{
      auto tmp = a*x_v(ss)+b*y_v(ss);
      coalescedWrite(inner_tmp_v[ss],innerProduct(tmp,tmp));
      coalescedWrite(z_v[ss],tmp);
  });

  nrm = real(TensorRemove(sumD_gpu(inner_tmp_v,sites)));
#else
  // CPU 
  typedef decltype(innerProductD(x_v[0],y_v[0])) inner_t;
  Vector<inner_t> inner_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
  
  accelerator_for( ss, sites, nsimd,{
      auto tmp = a*x_v(ss)+b*y_v(ss);
      inner_tmp_v[ss]=innerProductD(tmp,tmp);
      z_v[ss]=tmp;
  });
  // Already promoted to double
  nrm = real(TensorRemove(sum(inner_tmp_v,sites)));
#endif
  grid->GlobalSum(nrm);
  return nrm; 
}
 
template<class vobj> strong_inline void
innerProductNorm(ComplexD& ip, RealD &nrm, const Lattice<vobj> &left,const Lattice<vobj> &right)
{
  conformable(left,right);

  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_typeD vector_type;
  Vector<ComplexD> tmp(2);

  GridBase *grid = left.Grid();

  auto left_v=left.View();
  auto right_v=right.View();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();

#ifdef GRID_NVCC
  // GPU
  typedef decltype(innerProduct(left_v[0],right_v[0])) inner_t;
  typedef decltype(innerProduct(left_v[0],left_v[0])) norm_t;
  Vector<inner_t> inner_tmp(sites);
  Vector<norm_t> norm_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
  auto norm_tmp_v = &norm_tmp[0];

  accelerator_for( ss, sites, nsimd,{
      auto left_tmp = left_v(ss);
      coalescedWrite(inner_tmp_v[ss],innerProduct(left_tmp,right_v(ss)));
      coalescedWrite(norm_tmp_v[ss],innerProduct(left_tmp,left_tmp)));
  });

  tmp[0] = TensorRemove(sumD_gpu(inner_tmp_v,sites));
  tmp[1] = TensorRemove(sumD_gpu(norm_tmp_v,sites));
#else
  // CPU
  typedef decltype(innerProductD(left_v[0],right_v[0])) inner_t;
  typedef decltype(innerProductD(left_v[0],left_v[0])) norm_t;
  Vector<inner_t> inner_tmp(sites);
  Vector<norm_t> norm_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
  auto norm_tmp_v = &norm_tmp[0];

  accelerator_for( ss, sites, nsimd,{
      auto left_tmp = left_v(ss);
      inner_tmp_v[ss] = innerProductD(left_tmp,right_v(ss));
      norm_tmp_v[ss] = innerProductD(left_tmp,left_tmp);
  });
  // Already promoted to double
  tmp[0] = TensorRemove(sum(inner_tmp_v,sites));
  tmp[1] = TensorRemove(sum(norm_tmp_v,sites));
#endif
  grid->GlobalSumVector(&tmp[0],2); // keep norm Complex -> can use GlobalSumVector
  ip = tmp[0];
  nrm = real(tmp[1]);
}

template<class Op,class T1>
inline auto sum(const LatticeUnaryExpression<Op,T1> & expr)
  ->typename decltype(expr.op.func(eval(0,expr.arg1)))::scalar_object
{
  return sum(closure(expr));
}

template<class Op,class T1,class T2>
inline auto sum(const LatticeBinaryExpression<Op,T1,T2> & expr)
      ->typename decltype(expr.op.func(eval(0,expr.arg1),eval(0,expr.arg2)))::scalar_object
{
  return sum(closure(expr));
}


template<class Op,class T1,class T2,class T3>
inline auto sum(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr)
  ->typename decltype(expr.op.func(eval(0,expr.arg1),
				      eval(0,expr.arg2),
				      eval(0,expr.arg3)
				      ))::scalar_object
{
  return sum(closure(expr));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// sliceSum, sliceInnerProduct, sliceAxpy, sliceNorm etc...
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class vobj> inline void sliceSum(const Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
{
  ///////////////////////////////////////////////////////
  // FIXME precision promoted summation
  // may be important for correlation functions
  // But easily avoided by using double precision fields
  ///////////////////////////////////////////////////////
  typedef typename vobj::scalar_object sobj;
  GridBase  *grid = Data.Grid();
  assert(grid!=NULL);

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  assert(orthogdim >= 0);
  assert(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  Vector<vobj> lvSum(rd); // will locally sum vectors first
  Vector<sobj> lsSum(ld,Zero());                    // sum across these down to scalars
  ExtractBuffer<sobj> extracted(Nsimd);                  // splitting the SIMD

  result.resize(fd); // And then global sum to return the same vector to every node 
  for(int r=0;r<rd;r++){
    lvSum[r]=Zero();
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  // sum over reduced dimension planes, breaking out orthog dir
  // Parallel over orthog direction
  auto Data_v=Data.View();
  thread_for( r,rd, {
    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;
	lvSum[r]=lvSum[r]+Data_v[ss];
      }
    }
  });

  // Sum across simd lanes in the plane, breaking out orthog dir.
  Coordinate icoor(Nd);

  for(int rt=0;rt<rd;rt++){

    extract(lvSum[rt],extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx];

    }
  }
  
  // sum over nodes.
  sobj gsum;
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      gsum=lsSum[lt];
    } else {
      gsum=Zero();
    }

    grid->GlobalSum(gsum);

    result[t]=gsum;
  }
}

template<class vobj>
static void sliceInnerProductVector( std::vector<ComplexD> & result, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int orthogdim) 
{
  typedef typename vobj::vector_type   vector_type;
  typedef typename vobj::scalar_type   scalar_type;
  GridBase  *grid = lhs.Grid();
  assert(grid!=NULL);
  conformable(grid,rhs.Grid());

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  assert(orthogdim >= 0);
  assert(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  Vector<vector_type> lvSum(rd); // will locally sum vectors first
  Vector<scalar_type > lsSum(ld,scalar_type(0.0));                    // sum across these down to scalars
  ExtractBuffer<iScalar<scalar_type> > extracted(Nsimd);   // splitting the SIMD  

  result.resize(fd); // And then global sum to return the same vector to every node for IO to file
  for(int r=0;r<rd;r++){
    lvSum[r]=Zero();
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  auto lhv=lhs.View();
  auto rhv=rhs.View();
  thread_for( r,rd,{

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;
	vector_type vv = TensorRemove(innerProduct(lhv[ss],rhv[ss]));
	lvSum[r]=lvSum[r]+vv;
      }
    }
  });

  // Sum across simd lanes in the plane, breaking out orthog dir.
  Coordinate icoor(Nd);
  for(int rt=0;rt<rd;rt++){

    iScalar<vector_type> temp; 
    temp._internal = lvSum[rt];
    extract(temp,extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx]._internal;

    }
  }
  
  // sum over nodes.
  scalar_type gsum;
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      gsum=lsSum[lt];
    } else {
      gsum=scalar_type(0.0);
    }

    grid->GlobalSum(gsum);

    result[t]=gsum;
  }
}
template<class vobj>
static void sliceNorm (std::vector<RealD> &sn,const Lattice<vobj> &rhs,int Orthog) 
{
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  int Nblock = rhs.Grid()->GlobalDimensions()[Orthog];
  std::vector<ComplexD> ip(Nblock);
  sn.resize(Nblock);
  
  sliceInnerProductVector(ip,rhs,rhs,Orthog);
  for(int ss=0;ss<Nblock;ss++){
    sn[ss] = real(ip[ss]);
  }
};


template<class vobj>
static void sliceMaddVector(Lattice<vobj> &R,std::vector<RealD> &a,const Lattice<vobj> &X,const Lattice<vobj> &Y,
			    int orthogdim,RealD scale=1.0) 
{    
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::tensor_reduced tensor_reduced;
  
  scalar_type zscale(scale);

  GridBase *grid  = X.Grid();

  int Nsimd  =grid->Nsimd();
  int Nblock =grid->GlobalDimensions()[orthogdim];

  int fd     =grid->_fdimensions[orthogdim];
  int ld     =grid->_ldimensions[orthogdim];
  int rd     =grid->_rdimensions[orthogdim];

  int e1     =grid->_slice_nblock[orthogdim];
  int e2     =grid->_slice_block [orthogdim];
  int stride =grid->_slice_stride[orthogdim];

  Coordinate icoor;
  for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    vector_type    av;

    for(int l=0;l<Nsimd;l++){
      grid->iCoorFromIindex(icoor,l);
      int ldx =r+icoor[orthogdim]*rd;
      scalar_type *as =(scalar_type *)&av;
      as[l] = scalar_type(a[ldx])*zscale;
    }

    tensor_reduced at; at=av;

    auto Rv=R.View();
    auto Xv=X.View();
    auto Yv=Y.View();
    thread_for_collapse(2, n, e1, {
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;
	Rv[ss] = at*Xv[ss]+Yv[ss];
      }
    });
  }
};

/*
inline GridBase         *makeSubSliceGrid(const GridBase *BlockSolverGrid,int Orthog)
{
  int NN    = BlockSolverGrid->_ndimension;
  int nsimd = BlockSolverGrid->Nsimd();
  
  std::vector<int> latt_phys(0);
  std::vector<int> simd_phys(0);
  std::vector<int>  mpi_phys(0);
  
  for(int d=0;d<NN;d++){
    if( d!=Orthog ) { 
      latt_phys.push_back(BlockSolverGrid->_fdimensions[d]);
      simd_phys.push_back(BlockSolverGrid->_simd_layout[d]);
      mpi_phys.push_back(BlockSolverGrid->_processors[d]);
    }
  }
  return (GridBase *)new GridCartesian(latt_phys,simd_phys,mpi_phys); 
}
*/

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
  //  int nh =  FullGrid->_ndimension;
  //  int nl = SliceGrid->_ndimension;
  //  int nl = nh-1;

  //FIXME package in a convenient iterator
  //Should loop over a plane orthogonal to direction "Orthog"
  int stride=FullGrid->_slice_stride[Orthog];
  int block =FullGrid->_slice_block [Orthog];
  int nblock=FullGrid->_slice_nblock[Orthog];
  int ostride=FullGrid->_ostride[Orthog];

  auto X_v=X.View();
  auto Y_v=Y.View();
  auto R_v=R.View();
  thread_region
  {
    Vector<vobj> s_x(Nblock);

    thread_for_collapse_in_region(2, n,nblock, {
     for(int b=0;b<block;b++){
      int o  = n*stride + b;

      for(int i=0;i<Nblock;i++){
	s_x[i] = X_v[o+i*ostride];
      }

      vobj dot;
      for(int i=0;i<Nblock;i++){
	dot = Y_v[o+i*ostride];
	for(int j=0;j<Nblock;j++){
	  dot = dot + s_x[j]*(scale*aa(j,i));
	}
	R_v[o+i*ostride]=dot;
      }
    }});
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
  auto R_v = R.View();
  auto X_v = X.View();
  thread_region
  {
    std::vector<vobj> s_x(Nblock);


    thread_for_collapse_in_region( 2 ,n,nblock,{
    for(int b=0;b<block;b++){
      int o  = n*stride + b;

      for(int i=0;i<Nblock;i++){
	s_x[i] = X_v[o+i*ostride];
      }

      vobj dot;
      for(int i=0;i<Nblock;i++){
	dot = s_x[0]*(scale*aa(0,i));
	for(int j=1;j<Nblock;j++){
	  dot = dot + s_x[j]*(scale*aa(j,i));
	}
	R_v[o+i*ostride]=dot;
      }
    }});
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

  auto lhs_v=lhs.View();
  auto rhs_v=rhs.View();
  thread_region
  {
    std::vector<vobj> Left(Nblock);
    std::vector<vobj> Right(Nblock);
    Eigen::MatrixXcd  mat_thread = Eigen::MatrixXcd::Zero(Nblock,Nblock);

    thread_for_collapse_in_region( 2, n,nblock,{
    for(int b=0;b<block;b++){

      int o  = n*stride + b;

      for(int i=0;i<Nblock;i++){
	Left [i] = lhs_v[o+i*ostride];
	Right[i] = rhs_v[o+i*ostride];
      }

      for(int i=0;i<Nblock;i++){
      for(int j=0;j<Nblock;j++){
	auto tmp = innerProduct(Left[i],Right[j]);
	auto rtmp = TensorRemove(tmp);
	auto red  =  Reduce(rtmp);
	mat_thread(i,j) += std::complex<double>(real(red),imag(red));
      }}
    }});
    thread_critical
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




