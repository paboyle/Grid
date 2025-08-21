/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Cshift.h

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
#ifndef _GRID_FFT_H_
#define _GRID_FFT_H_

#ifdef GRID_CUDA
#include <cufft.h>
#endif

#ifdef GRID_HIP
#include <hipfft/hipfft.h>
#endif

#if !defined(GRID_CUDA) && !defined(GRID_HIP)
#ifdef HAVE_FFTW
#if defined(USE_MKL) || defined(GRID_SYCL)
#include <fftw/fftw3.h>
#else
#include <fftw3.h>
#endif
#endif
#endif

NAMESPACE_BEGIN(Grid);

#ifndef FFTW_FORWARD
#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)
#define FFTW_ESTIMATE (0)
#endif

template<class scalar> struct FFTW {
};

#ifdef GRID_HIP
template<> struct FFTW<ComplexD> {
public:
  static const int forward=FFTW_FORWARD;
  static const int backward=FFTW_BACKWARD;
  typedef hipfftDoubleComplex FFTW_scalar;
  typedef hipfftHandle        FFTW_plan;
  static FFTW_plan fftw_plan_many_dft(int rank, int *n,int howmany,
				      FFTW_scalar *in, int *inembed,		
				      int istride, int idist,		
				      FFTW_scalar *out, int *onembed,		
				      int ostride, int odist,		
				      int sign, unsigned flags) {
    FFTW_plan p;
    auto rv = hipfftPlanMany(&p,rank,n,n,istride,idist,n,ostride,odist,HIPFFT_Z2Z,howmany);
    GRID_ASSERT(rv==HIPFFT_SUCCESS);
    return p;
  }	  
    
  inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out, int sign) {
    hipfftResult rv;
    if ( sign == forward ) rv =hipfftExecZ2Z(p,in,out,HIPFFT_FORWARD);
    else                   rv =hipfftExecZ2Z(p,in,out,HIPFFT_BACKWARD);
    accelerator_barrier();
    GRID_ASSERT(rv==HIPFFT_SUCCESS);
  }
  inline static void fftw_destroy_plan(const FFTW_plan p) {
    hipfftDestroy(p);
  }
};
template<> struct FFTW<ComplexF> {
public:
  static const int forward=FFTW_FORWARD;
  static const int backward=FFTW_BACKWARD;
  typedef hipfftComplex      FFTW_scalar;
  typedef hipfftHandle        FFTW_plan;

  static FFTW_plan fftw_plan_many_dft(int rank, int *n,int howmany,
				      FFTW_scalar *in, int *inembed,		
				      int istride, int idist,		
				      FFTW_scalar *out, int *onembed,		
				      int ostride, int odist,		
				      int sign, unsigned flags) {
    FFTW_plan p;
    auto rv = hipfftPlanMany(&p,rank,n,n,istride,idist,n,ostride,odist,HIPFFT_C2C,howmany);
    GRID_ASSERT(rv==HIPFFT_SUCCESS);
    return p;
  }	  
    
  inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out, int sign) {
    hipfftResult rv;
    if ( sign == forward ) rv =hipfftExecC2C(p,in,out,HIPFFT_FORWARD);
    else                   rv =hipfftExecC2C(p,in,out,HIPFFT_BACKWARD);
    accelerator_barrier();
    GRID_ASSERT(rv==HIPFFT_SUCCESS);
  }
  inline static void fftw_destroy_plan(const FFTW_plan p) {
    hipfftDestroy(p);
  }
};
#endif

#ifdef GRID_CUDA
template<> struct FFTW<ComplexD> {
public:
  static const int forward=FFTW_FORWARD;
  static const int backward=FFTW_BACKWARD;
  typedef cufftDoubleComplex FFTW_scalar;
  typedef cufftHandle        FFTW_plan;

  static FFTW_plan fftw_plan_many_dft(int rank, int *n,int howmany,
				      FFTW_scalar *in, int *inembed,		
				      int istride, int idist,		
				      FFTW_scalar *out, int *onembed,		
				      int ostride, int odist,		
				      int sign, unsigned flags) {
    FFTW_plan p;
    cufftPlanMany(&p,rank,n,n,istride,idist,n,ostride,odist,CUFFT_Z2Z,howmany);
    return p;
  }	  
    
  inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out, int sign) {
    if ( sign == forward ) cufftExecZ2Z(p,in,out,CUFFT_FORWARD);
    else                   cufftExecZ2Z(p,in,out,CUFFT_BACKWARD);
    accelerator_barrier();
  }
  inline static void fftw_destroy_plan(const FFTW_plan p) {
    cufftDestroy(p);
  }
};
template<> struct FFTW<ComplexF> {
public:
  static const int forward=FFTW_FORWARD;
  static const int backward=FFTW_BACKWARD;
  typedef cufftComplex FFTW_scalar;
  typedef cufftHandle        FFTW_plan;

  static FFTW_plan fftw_plan_many_dft(int rank, int *n,int howmany,
				      FFTW_scalar *in, int *inembed,		
				      int istride, int idist,		
				      FFTW_scalar *out, int *onembed,		
				      int ostride, int odist,		
				      int sign, unsigned flags) {
    FFTW_plan p;
    cufftPlanMany(&p,rank,n,n,istride,idist,n,ostride,odist,CUFFT_C2C,howmany);
    return p;
  }	  
    
  inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out, int sign) {
    if ( sign == forward ) cufftExecC2C(p,in,out,CUFFT_FORWARD);
    else                   cufftExecC2C(p,in,out,CUFFT_BACKWARD);
    accelerator_barrier();
  }
  inline static void fftw_destroy_plan(const FFTW_plan p) {
    cufftDestroy(p);
  }
};
#endif

#if !defined(GRID_CUDA) && !defined(GRID_HIP)
#ifdef HAVE_FFTW
template<> struct FFTW<ComplexD> {
public:
  typedef fftw_complex FFTW_scalar;
  typedef fftw_plan    FFTW_plan;
  static FFTW_plan fftw_plan_many_dft(int rank, int *n,int howmany,
				      FFTW_scalar *in, int *inembed,		
				      int istride, int idist,		
				      FFTW_scalar *out, int *onembed,		
				      int ostride, int odist,		
				      int sign, unsigned flags) {
    return ::fftw_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
  }	  
    
  inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out, int sign) {
    ::fftw_execute_dft(p,in,out);
  }
  inline static void fftw_destroy_plan(const FFTW_plan p) {
    ::fftw_destroy_plan(p);
  }
};
template<> struct FFTW<ComplexF> {
public:
  typedef fftwf_complex FFTW_scalar;
  typedef fftwf_plan    FFTW_plan;
  static FFTW_plan fftw_plan_many_dft(int rank, int *n,int howmany,
				      FFTW_scalar *in, int *inembed,		
				      int istride, int idist,		
				      FFTW_scalar *out, int *onembed,		
				      int ostride, int odist,		
				      int sign, unsigned flags) {
    return ::fftwf_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
  }	  
    
  inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out, int sign) {
    ::fftwf_execute_dft(p,in,out);
  }
  inline static void fftw_destroy_plan(const FFTW_plan p) {
    ::fftwf_destroy_plan(p);
  }
};
#endif
#endif

class FFT {
private:
    
  double flops;
  double flops_call;
  uint64_t usec;
    
public:
    
  static const int forward=FFTW_FORWARD;
  static const int backward=FFTW_BACKWARD;
    
  double Flops(void) {return flops;}
  double MFlops(void) {return flops/usec;}
  double USec(void)   {return (double)usec;}    

  FFT ( GridCartesian * grid ) 
  {
    flops=0;
    usec =0;
  };
    
  ~FFT ( void)  {
    //    delete sgrid;
  }
    
  template<class vobj>
  void FFT_dim_mask(Lattice<vobj> &result,const Lattice<vobj> &source,Coordinate mask,int sign){

    //    vgrid=result.Grid();
    //    conformable(result.Grid(),vgrid);
    //    conformable(source.Grid(),vgrid);
    const int Ndim = source.Grid()->Nd();
    Lattice<vobj> tmp = source;
    for(int d=0;d<Ndim;d++){
      if( mask[d] ) {
	FFT_dim(result,tmp,d,sign);
	tmp=result;
      }
    }
  }

  template<class vobj>
  void FFT_all_dim(Lattice<vobj> &result,const Lattice<vobj> &source,int sign){
    const int Ndim = source.Grid()->Nd();
    Coordinate mask(Ndim,1);
    FFT_dim_mask(result,source,mask,sign);
  }


  template<class vobj>
  void FFT_dim(Lattice<vobj> &result,const Lattice<vobj> &source,int dim, int sign){
    const int Ndim = source.Grid()->Nd();
    GridBase *grid = source.Grid();
    conformable(result.Grid(),source.Grid());

    int L = grid->_ldimensions[dim];
    int G = grid->_fdimensions[dim];
      
    Coordinate layout(Ndim,1);
    
    // Construct pencils
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_type   scalar;
    typedef typename vobj::scalar_type   scalar_type;
    typedef typename vobj::vector_type   vector_type;
      
    //std::cout << "CPU view" << std::endl;
    
    typedef typename FFTW<scalar>::FFTW_scalar FFTW_scalar;
    typedef typename FFTW<scalar>::FFTW_plan   FFTW_plan;
      
    int Ncomp = sizeof(sobj)/sizeof(scalar);
    int64_t Nlow  = 1;
    int64_t Nhigh = 1;

    for(int d=0;d<dim;d++){
      Nlow*=grid->_ldimensions[d];
    }
    for(int d=dim+1;d<Ndim;d++){
      Nhigh*=grid->_ldimensions[d];
    }
    int64_t Nperp=Nlow*Nhigh;
    
    deviceVector<scalar> pgbuf; // Layout is [perp][component][dim]
    pgbuf.resize(Nperp*Ncomp*G);
    scalar *pgbuf_v = &pgbuf[0];
      
    int rank = 1;  /* 1d transforms */
    int n[] = {G}; /* 1d transforms of length G */
    int howmany = Ncomp * Nperp;
    int odist,idist,istride,ostride;
    idist   = odist   = G;            /* Distance between consecutive FT's */
    istride = ostride = 1;            /* Distance between two elements in the same FT */
    int *inembed = n, *onembed = n;
      
    scalar div;
    if ( sign == backward ) div = 1.0/G;
    else if ( sign == forward ) div = 1.0;
    else GRID_ASSERT(0);

    double t_pencil=0;
    double t_fft   =0;
    double t_total =-usecond();
    //    std::cout << GridLogPerformance<<"Making FFTW plan" << std::endl;
    /*
     *
     */
    FFTW_plan p;
    {
      FFTW_scalar *in = (FFTW_scalar *)&pgbuf_v[0];
      FFTW_scalar *out= (FFTW_scalar *)&pgbuf_v[0];
      p = FFTW<scalar>::fftw_plan_many_dft(rank,n,howmany,
					   in,inembed,
					   istride,idist,
					   out,onembed,
					   ostride, odist,
					   sign,FFTW_ESTIMATE);
    }
      
    // Barrel shift and collect global pencil
    //    std::cout << GridLogPerformance<<"Making pencil" << std::endl;
    Coordinate lcoor(Ndim), gcoor(Ndim);
    double t_copy=0;
    double t_shift=0;
    t_pencil = -usecond();
    result = source;
    int pc = grid->_processor_coor[dim];

    const Coordinate ldims = grid->_ldimensions;
    const Coordinate rdims = grid->_rdimensions;
    const Coordinate sdims = grid->_simd_layout;

    Coordinate processors = grid->_processors;
    Coordinate pgdims(Ndim);
    pgdims[0] = G;
    for(int d=0, dd=1;d<Ndim;d++){
      if ( d!=dim ) pgdims[dd++] = ldims[d];
    }
    int64_t pgvol=1;
    for(int d=0;d<Ndim;d++) pgvol*=pgdims[d];
    
    const int Nsimd = vobj::Nsimd();
    for(int p=0;p<processors[dim];p++) {
      t_copy-=usecond();
      autoView(r_v,result,AcceleratorRead);
      accelerator_for(idx, grid->oSites(), vobj::Nsimd(), {
#ifdef GRID_SIMT
      {
	int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
      for(int lane=0;lane<Nsimd;lane++) {
#endif
	Coordinate icoor;
	Coordinate ocoor;
	Coordinate pgcoor;

	Lexicographic::CoorFromIndex(icoor,lane,sdims);
	Lexicographic::CoorFromIndex(ocoor,idx,rdims);

	pgcoor[0] = ocoor[dim] + icoor[dim]*rdims[dim] + ((pc+p)%processors[dim])*L;
	for(int d=0,dd=1;d<Ndim;d++){
	  if ( d!=dim ) {
	    pgcoor[dd] = ocoor[d] + icoor[d]*rdims[d];
	    dd++;
	  }
	}

	// Map coordinates in lattice layout to FFTW index
	int64_t pgidx;
	Lexicographic::IndexFromCoor(pgcoor,pgidx,pgdims);

	vector_type *from = (vector_type *)&r_v[idx];
	scalar_type stmp;
	for(int w=0;w<Ncomp;w++){
	  int64_t pg_idx = pgidx + w*pgvol;
	  stmp = getlane(from[w], lane);
	  pgbuf_v[pg_idx] = stmp;
	}
#ifdef GRID_SIMT
      }
#else
      }
#endif
      });

      t_copy+=usecond();
      if (p != processors[dim] - 1) {
	Lattice<vobj> temp(grid);
	t_shift-=usecond();
	temp = Cshift(result,dim,L); result = temp;
	t_shift+=usecond();
      }
    }
    t_pencil += usecond();
      
    FFTW_scalar *in = (FFTW_scalar *)pgbuf_v;
    FFTW_scalar *out= (FFTW_scalar *)pgbuf_v;
    t_fft = -usecond();
    FFTW<scalar>::fftw_execute_dft(p,in,out,sign);
    t_fft += usecond();
    
    // performance counting
    flops_call = 5.0*howmany*G*log2(G);
    usec = t_fft;
    flops= flops_call;

    result = Zero();
    
    double t_insert = -usecond();
    {
      autoView(r_v,result,AcceleratorWrite);
      accelerator_for(idx,grid->oSites(),Nsimd,{
#ifdef GRID_SIMT
      {
	int lane=acceleratorSIMTlane(Nsimd); // buffer lane
#else
      for(int lane=0;lane<Nsimd;lane++) {
#endif
	Coordinate icoor(Ndim);
	Coordinate ocoor(Ndim);
	Coordinate pgcoor(Ndim);

	Lexicographic::CoorFromIndex(icoor,lane,sdims);
	Lexicographic::CoorFromIndex(ocoor,idx,rdims);

	pgcoor[0] = ocoor[dim] + icoor[dim]*rdims[dim] + pc*L;
	for(int d=0,dd=1;d<Ndim;d++){
	  if ( d!=dim ) {
	    pgcoor[dd] = ocoor[d] + icoor[d]*rdims[d];
	    dd++;
	  }
	}
	// Map coordinates in lattice layout to FFTW index
	int64_t pgidx;
	Lexicographic::IndexFromCoor(pgcoor,pgidx,pgdims);

	vector_type *to = (vector_type *)&r_v[idx];
	scalar_type stmp;
	for(int w=0;w<Ncomp;w++){
	  int64_t pg_idx = pgidx + w*pgvol;
	  stmp = pgbuf_v[pg_idx];
	  putlane(to[w], stmp, lane);
	}
	
#ifdef GRID_SIMT
      }
#else
      }
#endif
      });
    }

    result = result*div;

    t_insert +=usecond();
    
    // destroying plan
    FFTW<scalar>::fftw_destroy_plan(p);

    t_total +=usecond();

    std::cout <<GridLogPerformance<< " FFT took   "<<t_total/1.0e6 <<" s" << std::endl;
    std::cout <<GridLogPerformance<< " FFT pencil "<<t_pencil/1.0e6 <<" s" << std::endl;
    std::cout <<GridLogPerformance<< "  of which copy "<<t_copy/1.0e6 <<" s" << std::endl;
    std::cout <<GridLogPerformance<< "  of which shift"<<t_shift/1.0e6 <<" s" << std::endl;
    std::cout <<GridLogPerformance<< " FFT kernels "<<t_fft/1.0e6 <<" s" << std::endl;
    std::cout <<GridLogPerformance<< " FFT insert  "<<t_insert/1.0e6 <<" s" << std::endl;
    
  }
};

NAMESPACE_END(Grid);

#endif
