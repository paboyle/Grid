
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

#ifdef HAVE_FFTW
#ifdef USE_MKL
#include <fftw/fftw3.h>
#else
#include <fftw3.h>
#endif
#endif


namespace Grid {

  template<class scalar> struct FFTW { };

#ifdef HAVE_FFTW	
  template<> struct FFTW<ComplexD> {
  public:

    typedef fftw_complex FFTW_scalar;
    typedef fftw_plan    FFTW_plan;

    static FFTW_plan fftw_plan_many_dft(int rank, const int *n,int howmany,
					FFTW_scalar *in, const int *inembed,		
					int istride, int idist,		
					FFTW_scalar *out, const int *onembed,		
					int ostride, int odist,		
					int sign, unsigned flags) {
      return ::fftw_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
    }	  
    
    static void fftw_flops(const FFTW_plan p,double *add, double *mul, double *fmas){
      ::fftw_flops(p,add,mul,fmas);
    }

    inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out) {
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

    static FFTW_plan fftw_plan_many_dft(int rank, const int *n,int howmany,
					FFTW_scalar *in, const int *inembed,		
					int istride, int idist,		
					FFTW_scalar *out, const int *onembed,		
					int ostride, int odist,		
					int sign, unsigned flags) {
      return ::fftwf_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
    }	  
    
    static void fftw_flops(const FFTW_plan p,double *add, double *mul, double *fmas){
      ::fftwf_flops(p,add,mul,fmas);
    }

    inline static void fftw_execute_dft(const FFTW_plan p,FFTW_scalar *in,FFTW_scalar *out) {
      ::fftwf_execute_dft(p,in,out);
    }
    inline static void fftw_destroy_plan(const FFTW_plan p) {
      ::fftwf_destroy_plan(p);
    }
  };

#endif

#ifndef FFTW_FORWARD
#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)
#endif

  class FFT {
  private:
    
    GridCartesian *vgrid;
    GridCartesian *sgrid;
    
    int Nd;
    double flops;
    double flops_call;
    uint64_t usec;
    
    std::vector<int> dimensions;
    std::vector<int> processors;
    std::vector<int> processor_coor;
    
  public:
    
    static const int forward=FFTW_FORWARD;
    static const int backward=FFTW_BACKWARD;
    
    double Flops(void) {return flops;}
    double MFlops(void) {return flops/usec;}
    double USec(void)   {return (double)usec;}    

    FFT ( GridCartesian * grid ) :
    vgrid(grid),
    Nd(grid->_ndimension),
    dimensions(grid->_fdimensions),
    processors(grid->_processors),
    processor_coor(grid->_processor_coor)
    {
      flops=0;
      usec =0;
      std::vector<int> layout(Nd,1);
      sgrid = new GridCartesian(dimensions,layout,processors);
    };
    
    ~FFT ( void)  {
      delete sgrid;
    }
    
    template<class vobj>
    void FFT_dim_mask(Lattice<vobj> &result,const Lattice<vobj> &source,std::vector<int> mask,int sign){

      conformable(result._grid,vgrid);
      conformable(source._grid,vgrid);
      Lattice<vobj> tmp(vgrid);
      tmp = source;
      for(int d=0;d<Nd;d++){
	if( mask[d] ) {
	  FFT_dim(result,tmp,d,sign);
	  tmp=result;
	}
      }
    }

    template<class vobj>
    void FFT_all_dim(Lattice<vobj> &result,const Lattice<vobj> &source,int sign){
      std::vector<int> mask(Nd,1);
      FFT_dim_mask(result,source,mask,sign);
    }


    template<class vobj>
    void FFT_dim(Lattice<vobj> &result,const Lattice<vobj> &source,int dim, int sign){
#ifndef HAVE_FFTW
      assert(0);
#else
      conformable(result._grid,vgrid);
      conformable(source._grid,vgrid);

      int L = vgrid->_ldimensions[dim];
      int G = vgrid->_fdimensions[dim];
      
      std::vector<int> layout(Nd,1);
      std::vector<int> pencil_gd(vgrid->_fdimensions);
      
      pencil_gd[dim] = G*processors[dim];
      
      // Pencil global vol LxLxGxLxL per node
      GridCartesian pencil_g(pencil_gd,layout,processors);
      
      // Construct pencils
      typedef typename vobj::scalar_object sobj;
      typedef typename sobj::scalar_type   scalar;
      
      Lattice<sobj> pgbuf(&pencil_g);
      

      typedef typename FFTW<scalar>::FFTW_scalar FFTW_scalar;
      typedef typename FFTW<scalar>::FFTW_plan   FFTW_plan;
      
      int Ncomp = sizeof(sobj)/sizeof(scalar);
      int Nlow  = 1;
      for(int d=0;d<dim;d++){
        Nlow*=vgrid->_ldimensions[d];
      }
      
      int rank = 1;  /* 1d transforms */
      int n[] = {G}; /* 1d transforms of length G */
      int howmany = Ncomp;
      int odist,idist,istride,ostride;
      idist   = odist   = 1;          /* Distance between consecutive FT's */
      istride = ostride = Ncomp*Nlow; /* distance between two elements in the same FT */
      int *inembed = n, *onembed = n;
      
      scalar div;
	  if ( sign == backward ) div = 1.0/G;
	  else if ( sign == forward ) div = 1.0;
	  else assert(0);
      
      FFTW_plan p;
      {
        FFTW_scalar *in = (FFTW_scalar *)&pgbuf._odata[0];
        FFTW_scalar *out= (FFTW_scalar *)&pgbuf._odata[0];
        p = FFTW<scalar>::fftw_plan_many_dft(rank,n,howmany,
                                             in,inembed,
                                             istride,idist,
                                             out,onembed,
                                             ostride, odist,
                                             sign,FFTW_ESTIMATE);
      }
      
      // Barrel shift and collect global pencil
      std::vector<int> lcoor(Nd), gcoor(Nd);
      result = source;
      int pc = processor_coor[dim];
      for(int p=0;p<processors[dim];p++) {
        PARALLEL_REGION
        {
          std::vector<int> cbuf(Nd);
          sobj s;
          
          PARALLEL_FOR_LOOP_INTERN
          for(int idx=0;idx<sgrid->lSites();idx++) {
            sgrid->LocalIndexToLocalCoor(idx,cbuf);
            peekLocalSite(s,result,cbuf);
	    cbuf[dim]+=((pc+p) % processors[dim])*L;
	    //            cbuf[dim]+=p*L;
            pokeLocalSite(s,pgbuf,cbuf);
          }
        }
        if (p != processors[dim] - 1)
        {
          result = Cshift(result,dim,L);
        }
      }
      
      // Loop over orthog coords
      int NN=pencil_g.lSites();
      GridStopWatch timer;
      timer.Start();
      PARALLEL_REGION
      {
        std::vector<int> cbuf(Nd);
        
        PARALLEL_FOR_LOOP_INTERN
        for(int idx=0;idx<NN;idx++) {
          pencil_g.LocalIndexToLocalCoor(idx, cbuf);
          if ( cbuf[dim] == 0 ) {  // restricts loop to plane at lcoor[dim]==0
            FFTW_scalar *in = (FFTW_scalar *)&pgbuf._odata[idx];
            FFTW_scalar *out= (FFTW_scalar *)&pgbuf._odata[idx];
            FFTW<scalar>::fftw_execute_dft(p,in,out);
          }
        }
      }
      timer.Stop();
      
      // performance counting
      double add,mul,fma;
      FFTW<scalar>::fftw_flops(p,&add,&mul,&fma);
      flops_call = add+mul+2.0*fma;
      usec += timer.useconds();
      flops+= flops_call*NN;
      
      // writing out result
      PARALLEL_REGION
      {
        std::vector<int> clbuf(Nd), cgbuf(Nd);
        sobj s;
        
        PARALLEL_FOR_LOOP_INTERN
        for(int idx=0;idx<sgrid->lSites();idx++) {
          sgrid->LocalIndexToLocalCoor(idx,clbuf);
          cgbuf = clbuf;
          cgbuf[dim] = clbuf[dim]+L*pc;
          peekLocalSite(s,pgbuf,cgbuf);
          pokeLocalSite(s,result,clbuf);
        }
      }
      result = result*div;
      
      // destroying plan
      FFTW<scalar>::fftw_destroy_plan(p);
#endif
    }
  };
}

#endif
