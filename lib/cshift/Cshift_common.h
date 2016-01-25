    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/cshift/Cshift_common.h

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
#ifndef _GRID_CSHIFT_COMMON_H_
#define _GRID_CSHIFT_COMMON_H_

namespace Grid {

template<class vobj>
class SimpleCompressor {
public:
  void Point(int) {};

  vobj operator() (const vobj &arg,int dimension,int plane,int osite,GridBase *grid) {
    return arg;
  }
};

///////////////////////////////////////////////////////////////////
// Gather for when there is no need to SIMD split with compression
///////////////////////////////////////////////////////////////////
template<class vobj,class cobj,class compressor> void 
Gather_plane_simple (const Lattice<vobj> &rhs,std::vector<cobj,alignedAllocator<cobj> > &buffer,int dimension,int plane,int cbmask,compressor &compress, int off=0)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }
  
  int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  
  int e1=rhs._grid->_slice_nblock[dimension];
  int e2=rhs._grid->_slice_block[dimension];

  if ( cbmask == 0x3 ) { 
PARALLEL_NESTED_LOOP2
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o  = n*rhs._grid->_slice_stride[dimension];
	int bo = n*rhs._grid->_slice_block[dimension];
	buffer[off+bo+b]=compress(rhs._odata[so+o+b],dimension,plane,so+o+b,rhs._grid);
      }
    }
  } else { 
     int bo=0;
     for(int n=0;n<e1;n++){
       for(int b=0;b<e2;b++){
	 int o  = n*rhs._grid->_slice_stride[dimension];
	 int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
	 if ( ocb &cbmask ) {
	   buffer[off+bo++]=compress(rhs._odata[so+o+b],dimension,plane,so+o+b,rhs._grid);
	 }
       }
     }
  }
}


///////////////////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split with compression
///////////////////////////////////////////////////////////////////
template<class cobj,class vobj,class compressor> void 
Gather_plane_extract(const Lattice<vobj> &rhs,std::vector<typename cobj::scalar_object *> pointers,int dimension,int plane,int cbmask,compressor &compress)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }

  int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 

  int e1=rhs._grid->_slice_nblock[dimension];
  int e2=rhs._grid->_slice_block[dimension];
  
  if ( cbmask ==0x3){
PARALLEL_NESTED_LOOP2
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int o=n*rhs._grid->_slice_stride[dimension];
	int offset = b+n*rhs._grid->_slice_block[dimension];

	cobj temp =compress(rhs._odata[so+o+b],dimension,plane,so+o+b,rhs._grid);
	extract<cobj>(temp,pointers,offset);

      }
    }
  } else { 

    assert(0); //Fixme think this is buggy
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o=n*rhs._grid->_slice_stride[dimension];
	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);
	int offset = b+n*rhs._grid->_slice_block[dimension];

	if ( ocb & cbmask ) {
	  cobj temp =compress(rhs._odata[so+o+b],dimension,plane,so+o+b,rhs._grid);
	  extract<cobj>(temp,pointers,offset);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////
// Gather for when there is no need to SIMD split
//////////////////////////////////////////////////////
template<class vobj> void Gather_plane_simple (const Lattice<vobj> &rhs,std::vector<vobj,alignedAllocator<vobj> > &buffer,             int dimension,int plane,int cbmask)
{
  SimpleCompressor<vobj> dontcompress;
  Gather_plane_simple (rhs,buffer,dimension,plane,cbmask,dontcompress);
}

//////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split
//////////////////////////////////////////////////////
template<class vobj> void Gather_plane_extract(const Lattice<vobj> &rhs,std::vector<typename vobj::scalar_object *> pointers,int dimension,int plane,int cbmask)
{
  SimpleCompressor<vobj> dontcompress;
  Gather_plane_extract<vobj,vobj,decltype(dontcompress)>(rhs,pointers,dimension,plane,cbmask,dontcompress);
}

//////////////////////////////////////////////////////
// Scatter for when there is no need to SIMD split
//////////////////////////////////////////////////////
template<class vobj> void Scatter_plane_simple (Lattice<vobj> &rhs,std::vector<vobj,alignedAllocator<vobj> > &buffer, int dimension,int plane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    
  int e1=rhs._grid->_slice_nblock[dimension];
  int e2=rhs._grid->_slice_block[dimension];
  
  if ( cbmask ==0x3 ) {
PARALLEL_NESTED_LOOP2
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o   =n*rhs._grid->_slice_stride[dimension];
	int bo  =n*rhs._grid->_slice_block[dimension];
	rhs._odata[so+o+b]=buffer[bo+b];
      }
    }
  } else { 
    int bo=0;
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o   =n*rhs._grid->_slice_stride[dimension];
	int bo  =n*rhs._grid->_slice_block[dimension];
	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
	if ( ocb & cbmask ) {
	  rhs._odata[so+o+b]=buffer[bo++];
	}
      }
    }
  }
}

//////////////////////////////////////////////////////
// Scatter for when there *is* need to SIMD split
//////////////////////////////////////////////////////
 template<class vobj,class cobj> void Scatter_plane_merge(Lattice<vobj> &rhs,std::vector<cobj *> pointers,int dimension,int plane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
    
  int e1=rhs._grid->_slice_nblock[dimension];
  int e2=rhs._grid->_slice_block[dimension];

  if(cbmask ==0x3 ) {
PARALLEL_NESTED_LOOP2
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o      = n*rhs._grid->_slice_stride[dimension];
	int offset = b+n*rhs._grid->_slice_block[dimension];
	merge(rhs._odata[so+o+b],pointers,offset);
      }
    }
  } else { 
    assert(0); // think this is buggy FIXME
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o      = n*rhs._grid->_slice_stride[dimension];
	int offset = b+n*rhs._grid->_slice_block[dimension];
	int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);
	if ( ocb&cbmask ) {
	  merge(rhs._odata[so+o+b],pointers,offset);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////
// local to node block strided copies
//////////////////////////////////////////////////////
template<class vobj> void Copy_plane(Lattice<vobj>& lhs,const Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs._grid->_ostride[dimension]; // base offset for start of plane 

  int e1=rhs._grid->_slice_nblock[dimension]; // clearly loop invariant for icpc
  int e2=rhs._grid->_slice_block[dimension];

  if(cbmask == 0x3 ){
PARALLEL_NESTED_LOOP2
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
 
        int o =n*rhs._grid->_slice_stride[dimension]+b;
  	//lhs._odata[lo+o]=rhs._odata[ro+o];
	vstream(lhs._odata[lo+o],rhs._odata[ro+o]);
      }
    }
  } else { 
PARALLEL_NESTED_LOOP2
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
 
        int o =n*rhs._grid->_slice_stride[dimension]+b;
        int ocb=1<<lhs._grid->CheckerBoardFromOindex(o);
        if ( ocb&cbmask ) {
  	//lhs._odata[lo+o]=rhs._odata[ro+o];
	  vstream(lhs._odata[lo+o],rhs._odata[ro+o]);
	}
      }
    }
  }
  
}

template<class vobj> void Copy_plane_permute(Lattice<vobj>& lhs,const Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask,int permute_type)
{
 
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs._grid->_ostride[dimension]; // base offset for start of plane 

  int e1=rhs._grid->_slice_nblock[dimension];
  int e2=rhs._grid->_slice_block [dimension];
PARALLEL_NESTED_LOOP2
  for(int n=0;n<e1;n++){
  for(int b=0;b<e2;b++){

      int o  =n*rhs._grid->_slice_stride[dimension];
      int ocb=1<<lhs._grid->CheckerBoardFromOindex(o+b);
      if ( ocb&cbmask ) {
	permute(lhs._odata[lo+o+b],rhs._odata[ro+o+b],permute_type);
      }

  }}
}

//////////////////////////////////////////////////////
// Local to node Cshift
//////////////////////////////////////////////////////
template<class vobj> void Cshift_local(Lattice<vobj>& ret,const Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,Even);
  sshift[1] = rhs._grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,Odd);

  if ( sshift[0] == sshift[1] ) {
    Cshift_local(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_local(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_local(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> Lattice<vobj> Cshift_local(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid = rhs._grid;
  int fd = grid->_fdimensions[dimension];
  int rd = grid->_rdimensions[dimension];
  int ld = grid->_ldimensions[dimension];
  int gd = grid->_gdimensions[dimension];

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = grid->CheckerBoardDestination(rhs.checkerboard,shift,dimension);
  // the permute type
  int permute_dim =grid->PermuteDim(dimension);
  int permute_type=grid->PermuteType(dimension);

  for(int x=0;x<rd;x++){       

    int o   = 0;
    int bo  = x * grid->_ostride[dimension];
    
    int cb= (cbmask==0x2)? Odd : Even;

    int sshift = grid->CheckerBoardShiftForCB(rhs.checkerboard,dimension,shift,cb);
    int sx     = (x+sshift)%rd;

    int permute_slice=0;
    if(permute_dim){
      int wrap = sshift/rd;
      int  num = sshift%rd;
      if ( x< rd-num ) permute_slice=wrap;
      else permute_slice = 1-wrap;
    }

    if ( permute_slice ) Copy_plane_permute(ret,rhs,dimension,x,sx,cbmask,permute_type);
    else                 Copy_plane(ret,rhs,dimension,x,sx,cbmask); 

  
  }
  return ret;
}
}
#endif
