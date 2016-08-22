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

#include <Grid/fftw/fftw3.h>

namespace Grid {

  
  class FFT { 
  private:

    GridCartesian *vgrid;
    GridCartesian *sgrid;

    int Nd;
    std::vector<int> dimensions;
    std::vector<int> processors;
    std::vector<int> processor_coor;

  public:

    static const int forward=FFTW_FORWARD;
    static const int backward=FFTW_BACKWARD;

    FFT ( GridCartesian * grid ) : 
      vgrid(grid),
      Nd(grid->_ndimension),
      dimensions(grid->_fdimensions),
      processors(grid->_processors),
      processor_coor(grid->_processor_coor)
    {
      std::vector<int> layout(Nd,1);
      sgrid = new GridCartesian(dimensions,layout,processors);
    };

    ~FFT ( void)  { 
      delete sgrid; 
    }
    
    template<class vobj>
    void FFT_dim(Lattice<vobj> &result,const Lattice<vobj> &source,int dim, int inverse){

      conformable(result._grid,vgrid);
      conformable(source._grid,vgrid);

      int L = vgrid->_ldimensions[dim];
      int G = vgrid->_fdimensions[dim];

      std::vector<int> layout(Nd,1);
      std::vector<int> pencil_gd(vgrid->_fdimensions);
      std::vector<int> pencil_ld(processors);

      pencil_gd[dim] = G*processors[dim];    
      pencil_ld[dim] = G*processors[dim];    

      // Pencil global vol LxLxGxLxL per node
      GridCartesian pencil_g(pencil_gd,layout,processors);
      GridCartesian pencil_l(pencil_ld,layout,processors);

      // Construct pencils
      typedef typename vobj::scalar_object sobj;
      Lattice<vobj> ssource(vgrid); ssource =source;
      Lattice<sobj> pgsource(&pencil_g);
      Lattice<sobj> pgresult(&pencil_g);
      Lattice<sobj> plsource(&pencil_l);
      Lattice<sobj> plresult(&pencil_l);

      {
	assert(sizeof(typename sobj::scalar_type)==sizeof(ComplexD));
	assert(sizeof(fftw_complex)==sizeof(ComplexD));
	assert(sizeof(fftw_complex)==sizeof(ComplexD));

	int Ncomp = sizeof(sobj)/sizeof(fftw_complex);

	int rank = 1;  /* not 2: we are computing 1d transforms */
	int n[] = {G}; /* 1d transforms of length G */
	int howmany = Ncomp;
	int odist,idist,istride,ostride;
	idist   = odist   = 1;
	istride = ostride = Ncomp; /* distance between two elements in the same column */
	int *inembed = n, *onembed = n;

	fftw_complex *in = (fftw_complex *)&plsource._odata[0];
	fftw_complex *out= (fftw_complex *)&plresult._odata[0];
	
	int sign = FFTW_FORWARD;
	if (inverse) sign = FFTW_BACKWARD;

#ifdef HAVE_FFTW
	fftw_plan p = fftw_plan_many_dft(rank,n,howmany,
					 in,inembed,
					 istride,idist,
					 out,onembed,
					 ostride, odist,
					 sign,FFTW_ESTIMATE);
#else 
	fftw_plan p ;
	assert(0);
#endif

	// Barrel shift and collect global pencil
	for(int p=0;p<processors[dim];p++) { 

	  for(int idx=0;idx<sgrid->lSites();idx++) { 

	    std::vector<int> lcoor(Nd);
    	    sgrid->LocalIndexToLocalCoor(idx,lcoor);

	    sobj s;

	    peekLocalSite(s,ssource,lcoor);

	    lcoor[dim]+=p*L;
	   
	    pokeLocalSite(s,pgsource,lcoor);
	  }

	  ssource = Cshift(ssource,dim,L);
	}
	
	// Loop over orthog coords
	for(int idx=0;idx<sgrid->lSites();idx++) { 

	  std::vector<int> pcoor(Nd,0);
	  std::vector<int> lcoor(Nd);
	  sgrid->LocalIndexToLocalCoor(idx,lcoor);

	  if ( lcoor[dim] == 0 ) {  // restricts loop to plane at lcoor[dim]==0
	  
	    // Project to local pencil array
	    for(int l=0;l<G;l++){
	      sobj s;
	      pcoor[dim]=l;
	      lcoor[dim]=l;
	      peekLocalSite(s,pgsource,lcoor);
	      pokeLocalSite(s,plsource,pcoor);
	    }

	    // FFT the pencil
#ifdef HAVE_FFTW
	    fftw_execute(p);
#endif

	    // Extract the result
	    for(int l=0;l<L;l++){
	      sobj s;
	      int p = processor_coor[dim];
	      lcoor[dim] = l;
	      pcoor[dim] = l+L*p;
	      peekLocalSite(s,plresult,pcoor);
	      pokeLocalSite(s,result,lcoor);
	    }

	  }
	}
	  
	fftw_destroy_plan(p);
      }
    }

  };


}

#endif
