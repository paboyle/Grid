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

  std::vector<sobj> sumarray(nthread);
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

  std::vector<sobj> sumarray(nthread);
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

  bool ok;
#ifdef GRID_SYCL
  //  uint64_t csum=0;
  //  uint64_t csum2=0;
  //  if ( FlightRecorder::LoggingMode != FlightRecorder::LoggingModeNone)
  //  {
  // Hack
  // Fast integer xor checksum. Can also be used in comms now.
  //    autoView(l_v,left,AcceleratorRead);
  //    Integer words = left.Grid()->oSites()*sizeof(vobj)/sizeof(uint64_t);
  //    uint64_t *base= (uint64_t *)&l_v[0];
  //    csum=svm_xor(base,words);
  //    ok = FlightRecorder::CsumLog(csum);
  //    if ( !ok ) {
  //      csum2=svm_xor(base,words);
  //      std::cerr<< " Bad CSUM " << std::hex<< csum << " recomputed as "<<csum2<<std::dec<<std::endl;
  //    } else {
  //      csum2=svm_xor(base,words);
  //      std::cerr<< " ok CSUM " << std::hex<< csum << " recomputed as "<<csum2<<std::dec<<std::endl;
  //    }
  //    GRID_ASSERT(ok);
  // }
#endif
  FlightRecorder::StepLog("rank inner product");
  ComplexD nrm = rankInnerProduct(left,right);
  //  ComplexD nrmck=nrm;
  RealD local = real(nrm);
  ok = FlightRecorder::NormLog(real(nrm));
  if ( !ok ) {
    ComplexD nrm2 = rankInnerProduct(left,right);
    RealD local2 = real(nrm2);
    std::cerr<< " Bad NORM " << local << " recomputed as "<<local2<<std::endl;
    GRID_ASSERT(ok);
  }
  FlightRecorder::StepLog("Start global sum");
  grid->GlobalSumP2P(nrm);
  //  grid->GlobalSum(nrm);
  FlightRecorder::StepLog("Finished global sum");
  //  std::cout << " norm "<< nrm << " p2p norm "<<nrmck<<std::endl;
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
  typedef decltype(innerProduct(x_v[0],y_v[0])) inner_t;
  deviceVector<inner_t> inner_tmp;
  inner_tmp.resize(sites);
  auto inner_tmp_v = &inner_tmp[0];

  accelerator_for( ss, sites, nsimd,{
      auto tmp = a*x_v(ss)+b*y_v(ss);
      coalescedWrite(inner_tmp_v[ss],innerProduct(tmp,tmp));
      coalescedWrite(z_v[ss],tmp);
  });
  bool ok;
  nrm = real(TensorRemove(sumD(inner_tmp_v,sites)));
  ok = FlightRecorder::NormLog(real(nrm));
  GRID_ASSERT(ok);
  RealD local = real(nrm);
  grid->GlobalSum(nrm);
  FlightRecorder::ReductionLog(local,real(nrm));
  return nrm; 
}
 
template<class vobj> strong_inline void
innerProductNorm(ComplexD& ip, RealD &nrm, const Lattice<vobj> &left,const Lattice<vobj> &right)
{
  conformable(left,right);

  typedef typename vobj::vector_typeD vector_type;
  std::vector<ComplexD> tmp(2);

  GridBase *grid = left.Grid();

  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();

  // GPU
  typedef decltype(innerProductD(vobj(),vobj())) inner_t;
  typedef decltype(innerProductD(vobj(),vobj())) norm_t;
  deviceVector<inner_t> inner_tmp(sites);
  deviceVector<norm_t>  norm_tmp(sites);
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

template<class vobj> inline void sliceSum(const Lattice<vobj> &Data,
					  std::vector<typename vobj::scalar_object> &result,
					  int orthogdim)
{
  ///////////////////////////////////////////////////////
  // FIXME precision promoted summation
  // may be important for correlation functions
  // But easily avoided by using double precision fields
  ///////////////////////////////////////////////////////
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_object::scalar_type scalar_type;
  GridBase  *grid = Data.Grid();
  GRID_ASSERT(grid!=NULL);

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  GRID_ASSERT(orthogdim >= 0);
  GRID_ASSERT(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  std::vector<vobj> lvSum(rd); // will locally sum vectors first
  std::vector<sobj> lsSum(ld,Zero());                    // sum across these down to scalars
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
  //  std::cout << GridLogMessage << " sliceSum local"<<t_sum<<" us, host+mpi "<<t_rest<<std::endl;
  
}
template<class vobj> inline
std::vector<typename vobj::scalar_object> 
sliceSum(const Lattice<vobj> &Data,int orthogdim)
{
  std::vector<typename vobj::scalar_object> result;
  sliceSum(Data,result,orthogdim);
  return result;
}

/*
Reimplement

1)
template<class vobj>
static void sliceMaddMatrix (Lattice<vobj> &R,Eigen::MatrixXcd &aa,const Lattice<vobj> &X,const Lattice<vobj> &Y,int Orthog,RealD scale=1.0) 

2)
template<class vobj>
static void sliceInnerProductMatrix(  Eigen::MatrixXcd &mat, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int Orthog) 

3)
-- Make Slice Mul Matrix call sliceMaddMatrix
 */
template<class vobj>
static void sliceInnerProductVector( std::vector<ComplexD> & result, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int orthogdim) 
{
  typedef typename vobj::vector_type   vector_type;
  typedef typename vobj::scalar_type   scalar_type;
  GridBase  *grid = lhs.Grid();
  GRID_ASSERT(grid!=NULL);
  conformable(grid,rhs.Grid());

  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  GRID_ASSERT(orthogdim >= 0);
  GRID_ASSERT(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  std::vector<vector_type> lvSum(rd); // will locally sum vectors first
  std::vector<scalar_type > lsSum(ld,scalar_type(0.0));                    // sum across these down to scalars
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

inline GridBase         *makeSubSliceGrid(const GridBase *BlockSolverGrid,int Orthog)
{
  int NN    = BlockSolverGrid->_ndimension;
  int nsimd = BlockSolverGrid->Nsimd();
  
  std::vector<int> latt_phys(NN-1);
  Coordinate simd_phys;
  std::vector<int>  mpi_phys(NN-1);
  Coordinate checker_dim_mask(NN-1);
  int checker_dim=-1;

  int dd;
  for(int d=0;d<NN;d++){
    if( d!=Orthog ) { 
      latt_phys[dd]=BlockSolverGrid->_fdimensions[d];
      mpi_phys[dd] =BlockSolverGrid->_processors[d];
      checker_dim_mask[dd] = BlockSolverGrid->_checker_dim_mask[d];
      if ( d == BlockSolverGrid->_checker_dim ) checker_dim = dd;
      dd++;
    }
  }
  simd_phys=GridDefaultSimd(latt_phys.size(),nsimd);
  GridCartesian *tmp         = new GridCartesian(latt_phys,simd_phys,mpi_phys);
  if(BlockSolverGrid->_isCheckerBoarded) {
    GridRedBlackCartesian *ret = new GridRedBlackCartesian(tmp,checker_dim_mask,checker_dim);
    delete tmp;
    return (GridBase *) ret;
  } else { 
    return (GridBase *) tmp;
  }
}

template<class vobj>
static void sliceMaddMatrix (Lattice<vobj> &R,Eigen::MatrixXcd &aa,const Lattice<vobj> &X,const Lattice<vobj> &Y,int Orthog,RealD scale=1.0) 
{    
  GridBase *FullGrid = X.Grid();
  GridBase *SliceGrid = makeSubSliceGrid(FullGrid,Orthog);

  Lattice<vobj> Ys(SliceGrid);
  Lattice<vobj> Rs(SliceGrid);
  Lattice<vobj> Xs(SliceGrid);
  Lattice<vobj> RR(FullGrid);

  RR = R; // Copies checkerboard for insert
  
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::vector_type vector_type;
  int Nslice = X.Grid()->GlobalDimensions()[Orthog];
  for(int i=0;i<Nslice;i++){
    ExtractSlice(Ys,Y,i,Orthog);
    ExtractSlice(Rs,R,i,Orthog);
    Rs=Ys;
    for(int j=0;j<Nslice;j++){
      ExtractSlice(Xs,X,j,Orthog);
      Rs = Rs + Xs*(scale*aa(j,i));
    }
    InsertSlice(Rs,RR,i,Orthog);
  }
  R=RR; // Copy back handles arguments aliasing case
  delete SliceGrid;
};

template<class vobj>
static void sliceMulMatrix (Lattice<vobj> &R,Eigen::MatrixXcd &aa,const Lattice<vobj> &X,int Orthog,RealD scale=1.0)
{
  R=Zero();
  sliceMaddMatrix(R,aa,X,R,Orthog,scale);
};


template<class vobj>
static void sliceInnerProductMatrix(  Eigen::MatrixXcd &mat, const Lattice<vobj> &lhs,const Lattice<vobj> &rhs,int Orthog) 
{
  GridBase *SliceGrid = makeSubSliceGrid(lhs.Grid(),Orthog);

  Lattice<vobj> ls(SliceGrid);
  Lattice<vobj> rs(SliceGrid);
  
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::vector_type vector_type;
  int Nslice = lhs.Grid()->GlobalDimensions()[Orthog];
  mat = Eigen::MatrixXcd::Zero(Nslice,Nslice);
  for(int s=0;s<Nslice;s++){
    ExtractSlice(ls,lhs,s,Orthog);
    for(int ss=0;ss<Nslice;ss++){
      ExtractSlice(rs,rhs,ss,Orthog);
      mat(s,ss) = innerProduct(ls,rs);
    }
  }
  delete SliceGrid;
}

/////////////////////////////////////////////////////////////////////////////
// Batched linear algebra for Krylov orthogonalisation (GCR history windows).
//
//   innerProductMulti(out,left,right):  out[j] = <left[j],right> for all j in
//        ONE kernel (right read once), ONE device reduction / sync (the per-site
//        partials are an iVector over the batch) and ONE GlobalSumVector.
//   axpyMulti(z,b,x):                   z = z + sum_j b[j] x[j] in ONE pass.
//   axpyMultiNorm(z,b,x):               same, returning global |z|^2 from the
//        same pass (one reduction, one GlobalSum).
//
// The batch width is a template parameter chosen at runtime from {2,4,8,16}
// so a short window (mmax=2 smoother) does not pay for 16 partial lanes;
// windows longer than 16 are processed in chunks of 16, one reduction each.
// Same code path for every Lattice<vobj>: fine fermion fields and coarse
// multi-RHS fields (nrhs folded into the grid) alike.
/////////////////////////////////////////////////////////////////////////////
// The batch of views and coefficients is passed to the kernel BY VALUE as
// lambda-captured kernel arguments (a ViewPack), so there is no per-call
// deviceVector allocation and no synchronous host->device memcpy of pointer
// tables: the only sync points are the reductions themselves.
// LatticeView has no default constructor, so the pack holds raw aligned
// storage and views are copied in bytewise (as basisRotateJ does through
// acceleratorPut); the struct is trivially copyable as a kernel argument.
template<class View,int B>
struct ViewPack {
  alignas(View) unsigned char raw[B*sizeof(View)];
  ComplexD b[B];
  accelerator_inline const View & v(int j) const { return reinterpret_cast<const View *>(raw)[j]; }
  void set(int j,const View &view){ memcpy(raw+j*sizeof(View),&view,sizeof(View)); }
};

template<int B,class vobj>
void rankInnerProductMultiChunk(ComplexD *out,int m,
				const std::vector<const Lattice<vobj>*> &left,
				const Lattice<vobj> &right)
{
  typedef decltype(innerProductD(vobj(),vobj())) inner_t;
  typedef iVector<inner_t,B> batch_t;
  typedef decltype(right.View(AcceleratorRead)) View;

  GRID_ASSERT(m>=1 && m<=B);
  GridBase *grid = right.Grid();
  const uint64_t sites = grid->oSites();

  std::vector<View> h_v; h_v.reserve(m);
  ViewPack<View,B> pack;
  for(int j=0;j<m;j++){
    conformable(*left[j],right);
    h_v.push_back(left[j]->View(AcceleratorRead));
    pack.set(j,h_v[j]);
  }
  for(int j=m;j<B;j++) pack.set(j,h_v[0]);   // valid but unused lanes

  deviceVector<batch_t> partial(sites);
  batch_t *partial_v = &partial[0];
  {
    autoView(right_v,right,AcceleratorRead);
    accelerator_for(ss,sites,1,{
	auto r = right_v[ss];
	batch_t acc;
	for(int j=0;j<B;j++){
	  if ( j<m ) acc._internal[j] = innerProductD(pack.v(j)[ss],r);
	  else       zeroit(acc._internal[j]);
	}
	partial_v[ss] = acc;
    });
  }
  for(int j=0;j<m;j++) h_v[j].ViewClose();

  auto res = sum(partial_v,sites);   // one reduction for the whole batch
  for(int j=0;j<m;j++) out[j] = TensorRemove(res._internal[j]);
}

template<class vobj>
void rankInnerProductMulti(std::vector<ComplexD> &out,
			   const std::vector<const Lattice<vobj>*> &left,
			   const Lattice<vobj> &right)
{
  int m = left.size();
  out.resize(m);
  for(int j0=0;j0<m;j0+=16){
    int mm = std::min(16,m-j0);
    std::vector<const Lattice<vobj>*> sub(left.begin()+j0,left.begin()+j0+mm);
    if      ( mm<=2 ) rankInnerProductMultiChunk<2> (&out[j0],mm,sub,right);
    else if ( mm<=4 ) rankInnerProductMultiChunk<4> (&out[j0],mm,sub,right);
    else if ( mm<=8 ) rankInnerProductMultiChunk<8> (&out[j0],mm,sub,right);
    else              rankInnerProductMultiChunk<16>(&out[j0],mm,sub,right);
  }
}

template<class vobj>
void innerProductMulti(std::vector<ComplexD> &out,
		       const std::vector<const Lattice<vobj>*> &left,
		       const Lattice<vobj> &right)
{
  rankInnerProductMulti(out,left,right);
  if ( out.size() ) right.Grid()->GlobalSumVector(&out[0],(int)out.size());
}

// z = z + sum_{j<m} b[j] x[j] for one chunk of at most B vectors; if do_norm,
// also writes per-site |z|^2 into inner_tmp_v.
template<int B,class vobj>
void axpyMultiChunk(Lattice<vobj> &z,const ComplexD *b,
		    const std::vector<const Lattice<vobj>*> &x,int m,
		    int do_norm,
		    decltype(innerProduct(vobj(),vobj())) *inner_tmp_v)
{
  typedef decltype(z.View(AcceleratorRead)) View;
  GRID_ASSERT(m>=1 && m<=B);
  GridBase *grid = z.Grid();
  const uint64_t nsimd = grid->Nsimd();
  const uint64_t sites = grid->oSites();

  std::vector<View> h_v; h_v.reserve(m);
  ViewPack<View,B> pack;
  for(int j=0;j<m;j++){
    conformable(*x[j],z);
    GRID_ASSERT(x[j]!=&z);            // window must not alias the accumulator
    h_v.push_back(x[j]->View(AcceleratorRead));
    pack.set(j,h_v[j]);
    pack.b[j] = b[j];
  }
  for(int j=m;j<B;j++){ pack.set(j,h_v[0]); pack.b[j] = ComplexD(0.0); }

  autoView(z_v,z,AcceleratorWrite);
  accelerator_for(ss,sites,nsimd,{
      auto acc = coalescedRead(z_v[ss]);
      for(int j=0;j<B;j++) if ( j<m ) acc = acc + pack.b[j]*coalescedRead(pack.v(j)[ss]);
      coalescedWrite(z_v[ss],acc);
      if ( do_norm ) coalescedWrite(inner_tmp_v[ss],innerProduct(acc,acc));
  });
  for(int j=0;j<m;j++) h_v[j].ViewClose();
}

// z = z + sum_j b[j] x[j];  if do_norm, returns global |z|^2 from the same
// (last) pass.  Chunks of 16 for windows longer than 16.
template<class vobj>
RealD axpyMultiNormImpl(Lattice<vobj> &z,const std::vector<ComplexD> &b,
			const std::vector<const Lattice<vobj>*> &x,int do_norm)
{
  typedef decltype(innerProduct(vobj(),vobj())) inner_t;
  int m = x.size();
  GRID_ASSERT((int)b.size()>=m);
  GridBase *grid = z.Grid();
  const uint64_t sites = grid->oSites();

  deviceVector<inner_t> inner_tmp(do_norm ? sites : 1);
  inner_t *inner_tmp_v = &inner_tmp[0];

  if ( m==0 ) {
    if ( do_norm ) return norm2(z);
    return 0.0;
  }
  for(int j0=0;j0<m;j0+=16){
    int mm   = std::min(16,m-j0);
    int last = (j0+mm>=m);
    std::vector<const Lattice<vobj>*> sub(x.begin()+j0,x.begin()+j0+mm);
    int dn = do_norm && last;
    if      ( mm<=2 ) axpyMultiChunk<2> (z,&b[j0],sub,mm,dn,inner_tmp_v);
    else if ( mm<=4 ) axpyMultiChunk<4> (z,&b[j0],sub,mm,dn,inner_tmp_v);
    else if ( mm<=8 ) axpyMultiChunk<8> (z,&b[j0],sub,mm,dn,inner_tmp_v);
    else              axpyMultiChunk<16>(z,&b[j0],sub,mm,dn,inner_tmp_v);
  }

  RealD nrm = 0.0;
  if ( do_norm ) {
    nrm = real(TensorRemove(sumD(inner_tmp_v,sites)));
    grid->GlobalSum(nrm);
  }
  return nrm;
}
template<class vobj>
void axpyMulti(Lattice<vobj> &z,const std::vector<ComplexD> &b,const std::vector<const Lattice<vobj>*> &x)
{
  axpyMultiNormImpl(z,b,x,0);
}
template<class vobj>
RealD axpyMultiNorm(Lattice<vobj> &z,const std::vector<ComplexD> &b,const std::vector<const Lattice<vobj>*> &x)
{
  return axpyMultiNormImpl(z,b,x,1);
}

NAMESPACE_END(Grid);




