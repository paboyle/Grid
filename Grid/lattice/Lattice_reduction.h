/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 
    Source file: ./lib/lattice/Lattice_reduction.h
    Copyright (C) 2015
Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>
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


#if defined(GRID_CUDA)||defined(GRID_HIP)
#include <Grid/lattice/Lattice_reduction_gpu.h>
#endif
#if defined(GRID_SYCL)
#include <Grid/lattice/Lattice_reduction_sycl.h>
#endif
#include <Grid/lattice/Lattice_slicesum_core.h>

NAMESPACE_BEGIN(Grid);

//////////////////////////////////////////////////////
// FIXME this should promote to double and accumulate
//////////////////////////////////////////////////////
template<class vobj>
inline typename vobj::scalar_object sum_cpu(const vobj *arg, Integer osites)
{
  typedef typename vobj::scalar_object  sobj;

  //  const int Nsimd = vobj::Nsimd();
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
inline typename vobj::scalar_objectD sumD_cpu(const vobj *arg, Integer osites)
{
  typedef typename vobj::scalar_objectD  sobj;

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
/*
Threaded max, don't use for now
template<class Double>
inline Double max(const Double *arg, Integer osites)
{
  //  const int Nsimd = vobj::Nsimd();
  const int nthread = GridThread::GetThreads();

  std::vector<Double> maxarray(nthread);
  
  thread_for(thr,nthread, {
    int nwork, mywork, myoff;
    nwork = osites;
    GridThread::GetWork(nwork,thr,mywork,myoff);
    Double max=arg[0];
    for(int ss=myoff;ss<mywork+myoff; ss++){
      if( arg[ss] > max ) max = arg[ss];
    }
    maxarray[thr]=max;
  });
  
  Double tmax=maxarray[0];
  for(int i=0;i<nthread;i++){
    if (maxarray[i]>tmax) tmax = maxarray[i];
  } 
  return tmax;
}
*/
template<class vobj>
inline typename vobj::scalar_object sum(const vobj *arg, Integer osites)
{
#if defined(GRID_CUDA)||defined(GRID_HIP)||defined(GRID_SYCL)
  return sum_gpu(arg,osites);
#else
  return sum_cpu(arg,osites);
#endif  
}
template<class vobj>
inline typename vobj::scalar_objectD sumD(const vobj *arg, Integer osites)
{
#if defined(GRID_CUDA)||defined(GRID_HIP)||defined(GRID_SYCL)
  return sumD_gpu(arg,osites);
#else
  return sumD_cpu(arg,osites);
#endif  
}
template<class vobj>
inline typename vobj::scalar_objectD sumD_large(const vobj *arg, Integer osites)
{
#if defined(GRID_CUDA)||defined(GRID_HIP)||defined(GRID_SYCL)
  return sumD_gpu_large(arg,osites);
#else
  return sumD_cpu(arg,osites);
#endif  
}

template<class vobj>
inline typename vobj::scalar_object rankSum(const Lattice<vobj> &arg)
{
  Integer osites = arg.Grid()->oSites();
#if defined(GRID_CUDA)||defined(GRID_HIP)||defined(GRID_SYCL)
  autoView( arg_v, arg, AcceleratorRead);
  return sum_gpu(&arg_v[0],osites);
#else
  autoView(arg_v, arg, CpuRead);
  return sum_cpu(&arg_v[0],osites);
#endif  
}

template<class vobj>
inline typename vobj::scalar_object sum(const Lattice<vobj> &arg)
{
  auto ssum = rankSum(arg);
  arg.Grid()->GlobalSum(ssum);
  return ssum;
}

template<class vobj>
inline typename vobj::scalar_object rankSumLarge(const Lattice<vobj> &arg)
{
#if defined(GRID_CUDA)||defined(GRID_HIP)||defined(GRID_SYCL)
  autoView( arg_v, arg, AcceleratorRead);
  Integer osites = arg.Grid()->oSites();
  return sum_gpu_large(&arg_v[0],osites);
#else
  autoView(arg_v, arg, CpuRead);
  Integer osites = arg.Grid()->oSites();
  return sum_cpu(&arg_v[0],osites);
#endif
}

template<class vobj>
inline typename vobj::scalar_object sum_large(const Lattice<vobj> &arg)
{
  auto ssum = rankSumLarge(arg);
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


template<class Op,class T1>
inline auto norm2(const LatticeUnaryExpression<Op,T1> & expr)  ->RealD
{
  return norm2(closure(expr));
}

template<class Op,class T1,class T2>
inline auto norm2(const LatticeBinaryExpression<Op,T1,T2> & expr)      ->RealD
{
  return norm2(closure(expr));
}


template<class Op,class T1,class T2,class T3>
inline auto norm2(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr)      ->RealD
{
  return norm2(closure(expr));
}


//The global maximum of the site norm2
template<class vobj> inline RealD maxLocalNorm2(const Lattice<vobj> &arg)
{
  typedef typename vobj::tensor_reduced vscalar;  //iScalar<iScalar<.... <vPODtype> > >
  typedef typename vscalar::scalar_object  scalar;   //iScalar<iScalar<.... <PODtype> > >

  Lattice<vscalar> inner = localNorm2(arg);

  auto grid = arg.Grid();

  RealD max;
  for(int l=0;l<grid->lSites();l++){
    Coordinate coor;
    scalar val;
    RealD r;
    grid->LocalIndexToLocalCoor(l,coor);
    peekLocalSite(val,inner,coor);
    r=real(TensorRemove(val));
    if( (l==0) || (r>max)){
      max=r;
    }
  }
  grid->GlobalMax(max);
  return max;
}

// Double inner product
template<class vobj>
inline ComplexD rankInnerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right)
{
  typedef typename vobj::vector_typeD vector_type;
  ComplexD  nrm;
  
  GridBase *grid = left.Grid();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();
  
  // Might make all code paths go this way.
  typedef decltype(innerProduct(vobj(),vobj())) inner_t;
  deviceVector<inner_t> inner_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
    
  {
    autoView( left_v , left, AcceleratorRead);
    autoView( right_v,right, AcceleratorRead);

    // GPU - SIMT lane compliance...
    accelerator_for( ss, sites, nsimd,{
	auto x_l = left_v(ss);
	auto y_l = right_v(ss);
	coalescedWrite(inner_tmp_v[ss],innerProduct(x_l,y_l));
    });
  }
  // This is in single precision and fails some tests
  auto anrm = sumD(inner_tmp_v,sites);  
  nrm = anrm;
  return nrm;
}


template<class vobj>
inline ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right) {
  GridBase *grid = left.Grid();

#ifdef GRID_SYCL
  uint64_t csum=0;
  if ( FlightRecorder::LoggingMode != FlightRecorder::LoggingModeNone)
  {
    // Hack
    // Fast integer xor checksum. Can also be used in comms now.
    autoView(l_v,left,AcceleratorRead);
    Integer words = left.Grid()->oSites()*sizeof(vobj)/sizeof(uint64_t);
    uint64_t *base= (uint64_t *)&l_v[0];
    csum=svm_xor(base,words);
  }
  FlightRecorder::CsumLog(csum);
#endif
  ComplexD nrm = rankInnerProduct(left,right);
  RealD local = real(nrm);
  FlightRecorder::NormLog(real(nrm)); 
  grid->GlobalSum(nrm);
  FlightRecorder::ReductionLog(local,real(nrm)); 
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

  //  typedef typename vobj::vector_typeD vector_type;
  RealD  nrm;
  
  GridBase *grid = x.Grid();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();
  
  // GPU
  autoView( x_v, x, AcceleratorRead);
  autoView( y_v, y, AcceleratorRead);
  autoView( z_v, z, AcceleratorWrite);
#if 0
  typedef decltype(innerProductD(x_v[0],y_v[0])) inner_t;
  Vector<inner_t> inner_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];

  accelerator_for( ss, sites, nsimd,{
      auto tmp = a*x_v(ss)+b*y_v(ss);
      coalescedWrite(inner_tmp_v[ss],innerProductD(tmp,tmp));
      coalescedWrite(z_v[ss],tmp);
  });
  nrm = real(TensorRemove(sum(inner_tmp_v,sites)));
#else
  typedef decltype(innerProduct(x_v[0],y_v[0])) inner_t;
  deviceVector<inner_t> inner_tmp;
  inner_tmp.resize(sites);
  auto inner_tmp_v = &inner_tmp[0];

  accelerator_for( ss, sites, nsimd,{
      auto tmp = a*x_v(ss)+b*y_v(ss);
      coalescedWrite(inner_tmp_v[ss],innerProduct(tmp,tmp));
      coalescedWrite(z_v[ss],tmp);
  });
  nrm = real(TensorRemove(sumD(inner_tmp_v,sites)));
#endif
  grid->GlobalSum(nrm);
  return nrm; 
}
 
template<class vobj> strong_inline void
innerProductNorm(ComplexD& ip, RealD &nrm, const Lattice<vobj> &left,const Lattice<vobj> &right)
{
  conformable(left,right);

  typedef typename vobj::vector_typeD vector_type;
  Vector<ComplexD> tmp(2);

  GridBase *grid = left.Grid();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();

  // GPU
  typedef decltype(innerProductD(vobj(),vobj())) inner_t;
  typedef decltype(innerProductD(vobj(),vobj())) norm_t;
  Vector<inner_t> inner_tmp(sites);
  Vector<norm_t>  norm_tmp(sites);
  auto inner_tmp_v = &inner_tmp[0];
  auto norm_tmp_v = &norm_tmp[0];
  {
    autoView(left_v,left, AcceleratorRead);
    autoView(right_v,right,AcceleratorRead);
    accelerator_for( ss, sites, 1,{
	auto left_tmp = left_v[ss];
	inner_tmp_v[ss]=innerProductD(left_tmp,right_v[ss]);
        norm_tmp_v [ss]=innerProductD(left_tmp,left_tmp);
      });
  }

  tmp[0] = TensorRemove(sum(inner_tmp_v,sites));
  tmp[1] = TensorRemove(sum(norm_tmp_v,sites));

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
  typedef typename vobj::scalar_object::scalar_type scalar_type;
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
  int ostride=grid->_ostride[orthogdim];
  
  //Reduce Data down to lvSum
  sliceSumReduction(Data,lvSum,rd, e1,e2,stride,ostride,Nsimd);

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
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      result[t]=lsSum[lt];
    } else {
      result[t]=Zero();
    }

  }
  scalar_type * ptr = (scalar_type *) &result[0];
  int words = fd*sizeof(sobj)/sizeof(scalar_type);
  grid->GlobalSumVector(ptr, words);
}
template<class vobj> inline
std::vector<typename vobj::scalar_object> 
sliceSum(const Lattice<vobj> &Data,int orthogdim)
{
  std::vector<typename vobj::scalar_object> result;
  sliceSum(Data,result,orthogdim);
  return result;
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

  autoView( lhv, lhs, CpuRead);
  autoView( rhv, rhs, CpuRead);
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
  // perhaps easier to just promote A to a field and use regular madd
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
      av.putlane(scalar_type(a[ldx])*zscale,l);
    }

    tensor_reduced at; at=av;

    autoView( Rv, R, CpuWrite);
    autoView( Xv, X, CpuRead);
    autoView( Yv, Y, CpuRead);
    thread_for2d( n, e1, b,e2, {
	int ss= so+n*stride+b;
	Rv[ss] = at*Xv[ss]+Yv[ss];
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

  autoView( X_v, X, CpuRead);
  autoView( Y_v, Y, CpuRead);
  autoView( R_v, R, CpuWrite);
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
  // thread_for2d_in_region
  //Should loop over a plane orthogonal to direction "Orthog"
  int stride=FullGrid->_slice_stride[Orthog];
  int block =FullGrid->_slice_block [Orthog];
  int nblock=FullGrid->_slice_nblock[Orthog];
  int ostride=FullGrid->_ostride[Orthog];
  autoView( R_v, R, CpuWrite);
  autoView( X_v, X, CpuRead);
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

  autoView( lhs_v, lhs, CpuRead);
  autoView( rhs_v, rhs, CpuRead);
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




