    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/test_Grid_jacobi.cc

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
#include "Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class vobj>
class LinearOperator {
public:
  operator () (const Lattice<vobj>&src,Lattice<vobj> &result) {};
};

template<class vobj>
class LinearOperatorJacobi : public LinearOperator<vobj>
{
  CartesianStencil *Stencil;
  GridBase *_grid;
  std::vector<vobj,alignedAllocator<vobj> >  comm_buf;
  LinearOperatorJacobi(GridCartesian *grid)    
  {
    _grid = grid;
    int npoint=9;
    std::vector<int> directions(npoint);
    std::vector<int> displacements(npoint);
    for(int mu=0;mu<4;mu++){
      for(int mp=0;mp<2;mp++){
	int dir = 2*mu+mp;
	directions[dir]   = mu;
	displacements[dir]= -1+2*mp;
      }
    }
    directions[8] = 0;
    displacements[8] = 0;

    Stencil = new CartesianStencil(grid,npoint,0,directions,displacements);  

    comm_buf.resize(Stencil->_unified_buffer_size);

  }
  operator () (const Lattice<vobj>&src,Lattice<vobj> &result)
  {
    const int npoint=9;
    printf("calling halo exchange\n");fflush(stdout);
    myStencil.HaloExchange(Foo,comm_buf);
    vobj tmp;
    vobj
    for(int i=0;i<_grid->oSites();i++){
      for(int p=0;p<npoint;p++){

	int offset = Stencil->_offsets [p][i];
	int  local = Stencil->_is_local[p][i];
	int ptype  = Stencil->_permute_type[p];
	int perm   = Stencil->_permute[0][i];
	
	vobj *nbr;
	if ( local && perm ){
	  permute(tmp,src._odata[offset],ptype);
	  nbr = &tmp;
	} else if (local) {
	  nbr = &src._odata[offset];
	} else  {
	  nbr = &comm_buf[offset];
	}
       	result[i] = result[i]+*nbr;
      }
    }	
  }
  ~LinearOperatorJacobi()
  {
    delete Stencil;
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt_size  (4);
  std::vector<int> simd_layout(4);
  std::vector<int> mpi_layout (4);

  int omp=1;
  int lat=8;

  mpi_layout[0]=1;
  mpi_layout[1]=2;
  mpi_layout[2]=1;
  mpi_layout[3]=1;

    latt_size[0] = lat;
    latt_size[1] = lat;
    latt_size[2] = lat;
    latt_size[3] = lat;
    double volume = latt_size[0]*latt_size[1]*latt_size[2]*latt_size[3];
    
#ifdef AVX512
    simd_layout[0] = 1;
    simd_layout[1] = 2;
    simd_layout[2] = 2;
    simd_layout[3] = 2;
#endif
#if defined (AVX1)|| defined (AVX2)
    simd_layout[0] = 1;
    simd_layout[1] = 1;
    simd_layout[2] = 2;
    simd_layout[3] = 2;
#endif
#if defined (SSE2)
    simd_layout[0] = 1;
    simd_layout[1] = 1;
    simd_layout[2] = 1;
    simd_layout[3] = 2;
#endif
    
    GridCartesian Fine(latt_size,simd_layout,mpi_layout);
    GridRedBlackCartesian rbFine(latt_size,simd_layout,mpi_layout);
    
    LatticeColourMatrix Foo(&Fine);
    LatticeColourMatrix Bar(&Fine);
    LatticeColourMatrix Check(&Fine);
    LatticeColourMatrix Diff(&Fine);
    
    random(Foo);
    gaussian(Bar);


    for(int dir=0;dir<4;dir++){
      for(int disp=0;disp<Fine._rdimensions[dir];disp++){

	// start to test the Cartesian npoint stencil infrastructure
	int npoint=;
	std::vector<int> directions(npoint,dir);
	std::vector<int> displacements(npoint,disp);

	CartesianStencil myStencil(&Fine,npoint,0,directions,displacements);

	printf("STENCIL: osites %d %d dir %d disp %d\n",Fine.oSites(),(int)myStencil._offsets[0].size(),dir,disp);
	std::vector<int> ocoor(4);
	for(int o=0;o<Fine.oSites();o++){
	  Fine.oCoorFromOindex(ocoor,o);
	  ocoor[dir]=(ocoor[dir]+disp)%Fine._rdimensions[dir];
	  int nbr = Fine.oIndexReduced(ocoor);
	  int stcl= myStencil._offsets[0][o];
	  if(nbr!=stcl){
	    printf("STENCIL: nbr %d stencil._offset %d\n",nbr,stcl);
	  }
	}
	
	printf("allocating %d buffers\n",myStencil._unified_buffer_size);
	fflush(stdout);
	std::vector<vColourMatrix,alignedAllocator<vColourMatrix> >  comm_buf(myStencil._unified_buffer_size);
	printf("calling halo exchange\n");fflush(stdout);
	myStencil.HaloExchange(Foo,comm_buf);

	Bar = Cshift(Foo,dir,disp);

	// Implement a stencil code that should agree with cshift!
	for(int i=0;i<Check._grid->oSites();i++){

	  int offset = myStencil._offsets [0][i];
	  int  local = myStencil._is_local[0][i];
	  int permute_type = myStencil._permute_type[0];
	  int perm =myStencil._permute[0][i];
	  if ( local && perm )
	    permute(Check._odata[i],Foo._odata[offset],permute_type);
	  else if (local)
	    Check._odata[i] = Foo._odata[offset];
	  else 
	    Check._odata[i] = comm_buf[offset];
	  

	}


	std::vector<int> coor(4);
	for(coor[3]=0;coor[3]<latt_size[3]/mpi_layout[3];coor[3]++){
	for(coor[2]=0;coor[2]<latt_size[2]/mpi_layout[2];coor[2]++){
	for(coor[1]=0;coor[1]<latt_size[1]/mpi_layout[1];coor[1]++){
	for(coor[0]=0;coor[0]<latt_size[0]/mpi_layout[0];coor[0]++){

	  Complex diff;
	  ColourMatrix check,bar;
	  peekSite(check,Check,coor);
	  peekSite(bar,Bar,coor);

	  for(int r=0;r<3;r++){
	  for(int c=0;c<3;c++){
            diff =check._internal._internal[r][c]-bar._internal._internal[r][c];
            double nn=real(conjugate(diff)*diff);
            if ( nn > 0 ){
	      printf("Coor (%d %d %d %d) \t rc %d%d \t %le %le %le\n",
		     coor[0],coor[1],coor[2],coor[3],r,c,
		     nn,
		     real(check._internal._internal[r][c]),
		     real(bar._internal._internal[r][c])
		     );
	    }
	  }}
	 
	}}}}

      }
    }

 Grid_finalize();
}
