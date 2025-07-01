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

extern std::vector<std::pair<int,int> > Cshift_table; 
extern deviceVector<std::pair<int,int> > Cshift_table_device; 
extern std::vector<int> Cshift_vector;
extern deviceVector<int> Cshift_vector_device; 

// Copy Cshift map object (table or vector) to device
template<class vobj> 
inline void MapCshiftCopy(std::vector<vobj> &Cshift_obj, deviceVector<vobj> &Cshift_obj_device)
{
  // GPU version only
  uint64_t sz=Cshift_obj.size();
  if (Cshift_obj_device.size()!=sz )    {
    Cshift_obj_device.resize(sz);
  }
  acceleratorCopyToDevice((void *)&Cshift_obj[0],
			  (void *)&Cshift_obj_device[0],
			  sizeof(Cshift_obj[0])*sz);

}

// Copy Cshift map object (table or vector) to device and return pointer to device copy
template<class vobj> 
inline vobj *MapCshift(std::vector<vobj> &Cshift_obj, deviceVector<vobj> &Cshift_obj_device)
{
  MapCshiftCopy<vobj>(Cshift_obj, Cshift_obj_device);

  return &Cshift_obj_device[0];
}

// Calculate Cshift_vector
template<class vobj> 
void CalculateCshiftVector(Lattice<vobj> &ret, const Lattice<vobj> &rhs, int dimension, int cbmask)
{
  GridBase *grid = rhs.Grid();

  if ( !grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }
  
  int e1=grid->_slice_nblock[dimension]; // clearly loop invariant for icpc
  int e2=grid->_slice_block[dimension];
  int stride = grid->_slice_stride[dimension];

  if (Cshift_vector.size() < e1*e2) Cshift_vector.resize(e1*e2); // Let it grow to biggest 

  int ent = 0;
  if(cbmask == 0x3 ){
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
        int o =n*stride+b;
        Cshift_vector[ent++] = o;
      }
    }
  } else { 
    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
        int o =n*stride+b;
        int ocb=1<<ret.Grid()->CheckerBoardFromOindex(o);
        if ( ocb&cbmask ) {
          Cshift_vector[ent++] = o;
        }
      }
    }
  }
  
  if (ent < Cshift_vector.size()) Cshift_vector.resize(ent); // trim vector to actual size (relevant for checkerboarded dimensions)
}


///////////////////////////////////////////////////////////////////
// Gather for when there is no need to SIMD split 
///////////////////////////////////////////////////////////////////
template<class vobj> void 
Gather_plane_simple (const Lattice<vobj> &rhs,deviceVector<vobj> &buffer,int dimension,int plane,int cbmask, int off=0)
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
    auto table = MapCshift<std::pair<int,int> >(Cshift_table, Cshift_table_device);
    autoView(rhs_v , rhs, AcceleratorRead);
    accelerator_for(i,ent,vobj::Nsimd(),{
	coalescedWrite(buffer_p[table[i].first],coalescedRead(rhs_v[table[i].second]));
    });
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
    autoView(rhs_v , rhs, AcceleratorRead);
    accelerator_for(nn,e1*e2,1,{
	int n = nn%e1;
	int b = nn/e1;
	int o      =   n*n1;
	int offset = b+n*e2;
	
	vobj temp =rhs_v[so+o+b];
	extract<vobj>(temp,pointers,offset);
      });
  } else { 
    Coordinate rdim=rhs.Grid()->_rdimensions;
    Coordinate cdm =rhs.Grid()->_checker_dim_mask;
    std::cout << " Dense packed buffer WARNING " <<std::endl; // Does this get called twice once for each cb?
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
  }
}

//////////////////////////////////////////////////////
// Scatter for when there is no need to SIMD split
//////////////////////////////////////////////////////
template<class vobj> void Scatter_plane_simple (Lattice<vobj> &rhs,deviceVector<vobj> &buffer, int dimension,int plane,int cbmask)
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
    auto table = MapCshift<std::pair<int,int> >(Cshift_table, Cshift_table_device);
    autoView( rhs_v, rhs, AcceleratorWrite);
    accelerator_for(i,ent,vobj::Nsimd(),{
	coalescedWrite(rhs_v[table[i].first],coalescedRead(buffer_p[table[i].second]));
    });
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
    autoView( rhs_v , rhs, AcceleratorWrite);
    accelerator_for(nn,e1*e2,1,{
	int n = nn%e1;
	int b = nn/e1;
	int o      = n*_slice_stride;
	int offset = b+n*_slice_block;
	merge(rhs_v[so+o+b],pointers,offset);
      });
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

//////////////////////////////////////////////////////
// local to node block strided copies
//////////////////////////////////////////////////////

template<class vobj> void Copy_plane(Lattice<vobj>& lhs,const Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask)
{


  int ro  = rplane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs.Grid()->_ostride[dimension]; // base offset for start of plane 

  auto table = &Cshift_vector_device[0];

  autoView(rhs_v , rhs, AcceleratorRead);
  autoView(lhs_v , lhs, AcceleratorWrite);
  accelerator_for(i,Cshift_vector.size(),vobj::Nsimd(),{
    coalescedWrite(lhs_v[table[i]+lo],coalescedRead(rhs_v[table[i]+ro]));
  });
}

template<class vobj> void Copy_plane_permute(Lattice<vobj>& lhs,const Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask,int permute_type)
{

  int ro  = rplane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs.Grid()->_ostride[dimension]; // base offset for start of plane 

  auto table = &Cshift_vector_device[0];
  autoView( rhs_v, rhs, AcceleratorRead);
  autoView( lhs_v, lhs, AcceleratorWrite);
  accelerator_for(i,Cshift_vector.size(),1,{
    permute(lhs_v[table[i]+lo],rhs_v[table[i]+ro],permute_type);
    });
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

  int cb= (cbmask==0x2)? Odd : Even;
  int sshift = grid->CheckerBoardShiftForCB(rhs.Checkerboard(),dimension,shift,cb);

  // the permute type
  ret.Checkerboard() = grid->CheckerBoardDestination(rhs.Checkerboard(),shift,dimension);
  int permute_dim =grid->PermuteDim(dimension);
  int permute_type=grid->PermuteType(dimension);
  int permute_type_dist;

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
  int wrap = sshift/rd; wrap=wrap % ly;
  int  num = sshift%rd;

  // Calculate Cshift_vector - it's the same for all slices
  CalculateCshiftVector<vobj>(ret, rhs, dimension, cbmask);
  // Copy it to the device
  MapCshiftCopy<int>(Cshift_vector, Cshift_vector_device);

  for(int x=0;x<rd;x++){       

    int sx     = (x+sshift)%rd;
    
    int permute_slice=0;
    if(permute_dim){
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

