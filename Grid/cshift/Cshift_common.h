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
#pragma once

NAMESPACE_BEGIN(Grid);

extern Vector<std::pair<int,int> > Cshift_table; 

///////////////////////////////////////////////////////////////////
// Gather for when there is no need to SIMD split 
///////////////////////////////////////////////////////////////////
template<class vobj> void 
Gather_plane_simple (const Lattice<vobj> &rhs,cshiftVector<vobj> &buffer,int dimension,int plane,int cbmask, int off=0)
{
  int rd = rhs.Grid()->_rdimensions[dimension];

  if ( !rhs.Grid()->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }
  
  int so=plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
  int e1=rhs.Grid()->_slice_nblock[dimension];
  int e2=rhs.Grid()->_slice_block[dimension];
  int ent = 0;

  if(Cshift_table.size()<e1*e2) Cshift_table.resize(e1*e2); // Let it grow to biggest

  int stride=rhs.Grid()->_slice_stride[dimension];

  if ( cbmask == 0x3 ) { 
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o  = n*stride;
	int bo = n*e2;
	Cshift_table[ent++] = std::pair<int,int>(off+bo+b,so+o+b);
      }
    }
  } else { 
     int bo=0;
     for(int n=0;n<e1;n++){
       for(int b=0;b<e2;b++){
	 int o  = n*stride;
	 int ocb=1<<rhs.Grid()->CheckerBoardFromOindex(o+b);
	 if ( ocb &cbmask ) {
	   Cshift_table[ent++]=std::pair<int,int> (off+bo++,so+o+b);
	 }
       }
     }
  }
  {
    auto buffer_p = & buffer[0];
    auto table = &Cshift_table[0];
#ifdef ACCELERATOR_CSHIFT    
    autoView(rhs_v , rhs, AcceleratorRead);
    accelerator_for(i,ent,vobj::Nsimd(),{
	coalescedWrite(buffer_p[table[i].first],coalescedRead(rhs_v[table[i].second]));
    });
#else
    autoView(rhs_v , rhs, CpuRead);
    thread_for(i,ent,{
      buffer_p[table[i].first]=rhs_v[table[i].second];
    });
#endif
  }
}

///////////////////////////////////////////////////////////////////
// Gather for when there *is* need to SIMD split 
///////////////////////////////////////////////////////////////////
template<class vobj> void 
Gather_plane_extract(const Lattice<vobj> &rhs,
		     ExtractPointerArray<typename vobj::scalar_object> pointers,
		     int dimension,int plane,int cbmask)
{
  int rd = rhs.Grid()->_rdimensions[dimension];

  if ( !rhs.Grid()->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }

  int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 

  int e1=rhs.Grid()->_slice_nblock[dimension];
  int e2=rhs.Grid()->_slice_block[dimension];
  int n1=rhs.Grid()->_slice_stride[dimension];

  if ( cbmask ==0x3){
#ifdef ACCELERATOR_CSHIFT
    autoView(rhs_v , rhs, AcceleratorRead);
    accelerator_for(nn,e1*e2,1,{
	int n = nn%e1;
	int b = nn/e1;
	int o      =   n*n1;
	int offset = b+n*e2;
	
	vobj temp =rhs_v[so+o+b];
	extract<vobj>(temp,pointers,offset);
      });
#else
    autoView(rhs_v , rhs, CpuRead);
    thread_for2d(n,e1,b,e2,{
	int o      =   n*n1;
	int offset = b+n*e2;
	
	vobj temp =rhs_v[so+o+b];
	extract<vobj>(temp,pointers,offset);
      });
#endif
  } else { 
    Coordinate rdim=rhs.Grid()->_rdimensions;
    Coordinate cdm =rhs.Grid()->_checker_dim_mask;
    std::cout << " Dense packed buffer WARNING " <<std::endl; // Does this get called twice once for each cb?
#ifdef ACCELERATOR_CSHIFT    
    autoView(rhs_v , rhs, AcceleratorRead);
    accelerator_for(nn,e1*e2,1,{
	int n = nn%e1;
	int b = nn/e1;

	Coordinate coor;

	int o=n*n1;
	int oindex = o+b;

       	int cb = RedBlackCheckerBoardFromOindex(oindex, rdim, cdm);

	int ocb=1<<cb;
	int offset = b+n*e2;

	if ( ocb & cbmask ) {
	  vobj temp =rhs_v[so+o+b];
	  extract<vobj>(temp,pointers,offset);
	}
      });
#else
    autoView(rhs_v , rhs, CpuRead);
    thread_for2d(n,e1,b,e2,{

	Coordinate coor;

	int o=n*n1;
	int oindex = o+b;

       	int cb = RedBlackCheckerBoardFromOindex(oindex, rdim, cdm);

	int ocb=1<<cb;
	int offset = b+n*e2;

	if ( ocb & cbmask ) {
	  vobj temp =rhs_v[so+o+b];
	  extract<vobj>(temp,pointers,offset);
	}
      });
#endif
  }
}

//////////////////////////////////////////////////////
// Scatter for when there is no need to SIMD split
//////////////////////////////////////////////////////
template<class vobj> void Scatter_plane_simple (Lattice<vobj> &rhs,cshiftVector<vobj> &buffer, int dimension,int plane,int cbmask)
{
  int rd = rhs.Grid()->_rdimensions[dimension];

  if ( !rhs.Grid()->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
    
  int e1=rhs.Grid()->_slice_nblock[dimension];
  int e2=rhs.Grid()->_slice_block[dimension];
  int stride=rhs.Grid()->_slice_stride[dimension];

  if(Cshift_table.size()<e1*e2) Cshift_table.resize(e1*e2); // Let it grow to biggest

  int ent    =0;

  if ( cbmask ==0x3 ) {

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o   =n*rhs.Grid()->_slice_stride[dimension];
	int bo  =n*rhs.Grid()->_slice_block[dimension];
	Cshift_table[ent++] = std::pair<int,int>(so+o+b,bo+b);
      }
    }

  } else { 
    int bo=0;
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o   =n*rhs.Grid()->_slice_stride[dimension];
	int ocb=1<<rhs.Grid()->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
	if ( ocb & cbmask ) {
	  Cshift_table[ent++]=std::pair<int,int> (so+o+b,bo++);
	}
      }
    }
  }
  
  {
    auto buffer_p = & buffer[0];
    auto table = &Cshift_table[0];
#ifdef ACCELERATOR_CSHIFT    
    autoView( rhs_v, rhs, AcceleratorWrite);
    accelerator_for(i,ent,vobj::Nsimd(),{
	coalescedWrite(rhs_v[table[i].first],coalescedRead(buffer_p[table[i].second]));
    });
#else
    autoView( rhs_v, rhs, CpuWrite);
    thread_for(i,ent,{
      rhs_v[table[i].first]=buffer_p[table[i].second];
    });
#endif
  }
}

//////////////////////////////////////////////////////
// Scatter for when there *is* need to SIMD split
//////////////////////////////////////////////////////
template<class vobj> void Scatter_plane_merge(Lattice<vobj> &rhs,ExtractPointerArray<typename vobj::scalar_object> pointers,int dimension,int plane,int cbmask)
{
  int rd = rhs.Grid()->_rdimensions[dimension];

  if ( !rhs.Grid()->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
    
  int e1=rhs.Grid()->_slice_nblock[dimension];
  int e2=rhs.Grid()->_slice_block[dimension];

  if(cbmask ==0x3 ) {
    int _slice_stride = rhs.Grid()->_slice_stride[dimension];
    int _slice_block = rhs.Grid()->_slice_block[dimension];
#ifdef ACCELERATOR_CSHIFT    
    autoView( rhs_v , rhs, AcceleratorWrite);
    accelerator_for(nn,e1*e2,1,{
	int n = nn%e1;
	int b = nn/e1;
	int o      = n*_slice_stride;
	int offset = b+n*_slice_block;
	merge(rhs_v[so+o+b],pointers,offset);
      });
#else
    autoView( rhs_v , rhs, CpuWrite);
    thread_for2d(n,e1,b,e2,{
	int o      = n*_slice_stride;
	int offset = b+n*_slice_block;
	merge(rhs_v[so+o+b],pointers,offset);
    });
#endif
  } else { 

    // Case of SIMD split AND checker dim cannot currently be hit, except in 
    // Test_cshift_red_black code.
    std::cout << "Scatter_plane merge assert(0); think this is buggy FIXME "<< std::endl;// think this is buggy FIXME
    std::cout<<" Unthreaded warning -- buffer is not densely packed ??"<<std::endl;
    assert(0); // This will fail if hit on GPU
    autoView( rhs_v, rhs, CpuWrite);
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int o      = n*rhs.Grid()->_slice_stride[dimension];
	int offset = b+n*rhs.Grid()->_slice_block[dimension];
	int ocb=1<<rhs.Grid()->CheckerBoardFromOindex(o+b);
	if ( ocb&cbmask ) {
	  merge(rhs_v[so+o+b],pointers,offset);
	}
      }
    }
  }
}

#if (defined(GRID_CUDA) || defined(GRID_HIP)) && defined(ACCELERATOR_CSHIFT)

template <typename T>
T iDivUp(T a, T b) // Round a / b to nearest higher integer value
{ return (a % b != 0) ? (a / b + 1) : (a / b); }

template <typename T>
__global__ void populate_Cshift_table(T* vector, T lo, T ro, T e1, T e2, T stride)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= e1*e2) return;

    int n, b, o;

    n = idx / e2;
    b = idx % e2;
    o = n*stride + b;

    vector[2*idx + 0] = lo + o;
    vector[2*idx + 1] = ro + o;
}

#endif

//////////////////////////////////////////////////////
// local to node block strided copies
//////////////////////////////////////////////////////

template<class vobj> void Copy_plane(Lattice<vobj>& lhs,const Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask)
{
  int rd = rhs.Grid()->_rdimensions[dimension];

  if ( !rhs.Grid()->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int ro  = rplane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs.Grid()->_ostride[dimension]; // base offset for start of plane 

  int e1=rhs.Grid()->_slice_nblock[dimension]; // clearly loop invariant for icpc
  int e2=rhs.Grid()->_slice_block[dimension];
  int stride = rhs.Grid()->_slice_stride[dimension];

  if(Cshift_table.size()<e1*e2) Cshift_table.resize(e1*e2); // Let it grow to biggest

  int ent=0;

  if(cbmask == 0x3 ){
#if (defined(GRID_CUDA) || defined(GRID_HIP)) && defined(ACCELERATOR_CSHIFT)
    ent = e1*e2;
    dim3 blockSize(acceleratorThreads());
    dim3 gridSize(iDivUp((unsigned int)ent, blockSize.x));
    populate_Cshift_table<<<gridSize, blockSize>>>(&Cshift_table[0].first, lo, ro, e1, e2, stride);
    accelerator_barrier();
#else
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
        int o =n*stride+b;
	Cshift_table[ent++] = std::pair<int,int>(lo+o,ro+o);
      }
    }
#endif
  } else { 
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
        int o =n*stride+b;
        int ocb=1<<lhs.Grid()->CheckerBoardFromOindex(o);
        if ( ocb&cbmask ) {
	  Cshift_table[ent++] = std::pair<int,int>(lo+o,ro+o);
	}
      }
    }
  }

  {
    auto table = &Cshift_table[0];
#ifdef ACCELERATOR_CSHIFT    
    autoView(rhs_v , rhs, AcceleratorRead);
    autoView(lhs_v , lhs, AcceleratorWrite);
    accelerator_for(i,ent,vobj::Nsimd(),{
      coalescedWrite(lhs_v[table[i].first],coalescedRead(rhs_v[table[i].second]));
    });
#else
    autoView(rhs_v , rhs, CpuRead);
    autoView(lhs_v , lhs, CpuWrite);
    thread_for(i,ent,{
      lhs_v[table[i].first]=rhs_v[table[i].second];
    });
#endif
  }
}

template<class vobj> void Copy_plane_permute(Lattice<vobj>& lhs,const Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask,int permute_type)
{
  int rd = rhs.Grid()->_rdimensions[dimension];

  if ( !rhs.Grid()->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int ro  = rplane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs.Grid()->_ostride[dimension]; // base offset for start of plane 

  int e1=rhs.Grid()->_slice_nblock[dimension];
  int e2=rhs.Grid()->_slice_block [dimension];
  int stride = rhs.Grid()->_slice_stride[dimension];

  if(Cshift_table.size()<e1*e2) Cshift_table.resize(e1*e2); // Let it grow to biggest

  int ent=0;

  if ( cbmask == 0x3 ) {
#if (defined(GRID_CUDA) || defined(GRID_HIP)) && defined(ACCELERATOR_CSHIFT)
    ent = e1*e2;
    dim3 blockSize(acceleratorThreads());
    dim3 gridSize(iDivUp((unsigned int)ent, blockSize.x));
    populate_Cshift_table<<<gridSize, blockSize>>>(&Cshift_table[0].first, lo, ro, e1, e2, stride);
    accelerator_barrier();
#else
    for(int n=0;n<e1;n++){
    for(int b=0;b<e2;b++){
      int o  =n*stride;
      Cshift_table[ent++] = std::pair<int,int>(lo+o+b,ro+o+b);
    }}
#endif
  } else {
    for(int n=0;n<e1;n++){
    for(int b=0;b<e2;b++){
      int o  =n*stride;
      int ocb=1<<lhs.Grid()->CheckerBoardFromOindex(o+b);
      if ( ocb&cbmask ) Cshift_table[ent++] = std::pair<int,int>(lo+o+b,ro+o+b);
    }}
  }

  {
    auto table = &Cshift_table[0];
#ifdef ACCELERATOR_CSHIFT    
    autoView( rhs_v, rhs, AcceleratorRead);
    autoView( lhs_v, lhs, AcceleratorWrite);
    accelerator_for(i,ent,1,{
      permute(lhs_v[table[i].first],rhs_v[table[i].second],permute_type);
    });
#else
    autoView( rhs_v, rhs, CpuRead);
    autoView( lhs_v, lhs, CpuWrite);
    thread_for(i,ent,{
      permute(lhs_v[table[i].first],rhs_v[table[i].second],permute_type);
    });
#endif
  }
}

//////////////////////////////////////////////////////
// Local to node Cshift
//////////////////////////////////////////////////////
template<class vobj> void Cshift_local(Lattice<vobj>& ret,const Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,Even);
  sshift[1] = rhs.Grid()->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,Odd);

  if ( sshift[0] == sshift[1] ) {
    Cshift_local(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_local(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_local(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> void Cshift_local(Lattice<vobj> &ret,const Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid = rhs.Grid();
  int fd = grid->_fdimensions[dimension];
  int rd = grid->_rdimensions[dimension];
  int ld = grid->_ldimensions[dimension];
  int gd = grid->_gdimensions[dimension];
  int ly = grid->_simd_layout[dimension];

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  // the permute type
  ret.Checkerboard() = grid->CheckerBoardDestination(rhs.Checkerboard(),shift,dimension);
  int permute_dim =grid->PermuteDim(dimension);
  int permute_type=grid->PermuteType(dimension);
  int permute_type_dist;

  for(int x=0;x<rd;x++){       

    //    int o   = 0;
    int bo  = x * grid->_ostride[dimension];
    int cb= (cbmask==0x2)? Odd : Even;

    int sshift = grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);
    int sx     = (x+sshift)%rd;
    
    // wrap is whether sshift > rd.
    //  num is sshift mod rd.
    // 
    //  shift 7
    //
    //  XoXo YcYc 
    //  oXoX cYcY
    //  XoXo YcYc
    //  oXoX cYcY
    //
    //  sshift -- 
    //
    //  XX YY ; 3
    //  XX YY ; 0
    //  XX YY ; 3
    //  XX YY ; 0
    //
    int permute_slice=0;
    if(permute_dim){
      int wrap = sshift/rd; wrap=wrap % ly;
      int  num = sshift%rd;

      if ( x< rd-num ) permute_slice=wrap;
      else permute_slice = (wrap+1)%ly;

      if ( (ly>2) && (permute_slice) ) {
	assert(permute_type & RotateBit);
	permute_type_dist = permute_type|permute_slice;
      } else {
	permute_type_dist = permute_type;
      }
    }

    if ( permute_slice ) Copy_plane_permute(ret,rhs,dimension,x,sx,cbmask,permute_type_dist);
    else                 Copy_plane(ret,rhs,dimension,x,sx,cbmask); 
  
  }
}
NAMESPACE_END(Grid);

