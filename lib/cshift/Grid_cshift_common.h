#ifndef _GRID_CSHIFT_COMMON_H_
#define _GRID_CSHIFT_COMMON_H_

namespace Grid {

template<class vobj>
class SimpleCompressor {
public:
  void Point(int) {};

  vobj operator() (const vobj &arg) {
    return arg;
  }
};

///////////////////////////////////////////////////////////////////
// Gather for when there is no need to SIMD split with compression
///////////////////////////////////////////////////////////////////
template<class vobj,class cobj,class compressor> void 
Gather_plane_simple (const Lattice<vobj> &rhs,std::vector<cobj,alignedAllocator<cobj> > &buffer,int dimension,int plane,int cbmask,compressor &compress)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask = 0x3;
  }

  int so  = plane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  int o   = 0;                                      // relative offset to base within plane
  int bo  = 0;                                      // offset in buffer
  
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
      int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
      if ( ocb &cbmask ) {
	buffer[bo]=compress(rhs._odata[so+o+b]);
	bo++;
      }
    }
    o +=rhs._grid->_slice_stride[dimension];
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
  int o   = 0;                                      // relative offset to base within plane
  int bo  = 0;                                      // offset in buffer
    
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

      int offset = b+n*rhs._grid->_slice_block[dimension];

      int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);
      if ( ocb & cbmask ) {
	cobj temp; 
	temp =compress(rhs._odata[so+o+b]);
	extract<cobj>(temp,pointers,offset);
      }
      
    }
    o +=rhs._grid->_slice_stride[dimension];
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
  int o   = 0;                                      // relative offset to base within plane
  int bo  = 0;                                      // offset in buffer
    
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

      int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);// Could easily be a table lookup
      if ( ocb & cbmask ) {
	rhs._odata[so+o+b]=buffer[bo++];
      }

    }
    o +=rhs._grid->_slice_stride[dimension];
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
  int o   = 0;                                      // relative offset to base within plane
  int bo  = 0;                                      // offset in buffer
    
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){

      int offset = b+n*rhs._grid->_slice_block[dimension];
      int ocb=1<<rhs._grid->CheckerBoardFromOindex(o+b);
      if ( ocb&cbmask ) {
	merge(rhs._odata[so+o+b],pointers,offset);
      }
      
    }
    o +=rhs._grid->_slice_stride[dimension];
  }
}

//////////////////////////////////////////////////////
// local to node block strided copies
//////////////////////////////////////////////////////
template<class vobj> void Copy_plane(Lattice<vobj>& lhs,Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }


  int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs._grid->_ostride[dimension]; // base offset for start of plane 
  int o   = 0;                                     // relative offset to base within plane
  
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
      
      int ocb=1<<lhs._grid->CheckerBoardFromOindex(o+b);
      
      if ( ocb&cbmask ) {
	lhs._odata[lo+o+b]=rhs._odata[ro+o+b];
      }
      
    }
    o +=rhs._grid->_slice_stride[dimension];
  }
}

template<class vobj> void Copy_plane_permute(Lattice<vobj>& lhs,Lattice<vobj> &rhs, int dimension,int lplane,int rplane,int cbmask,int permute_type)
{
  int rd = rhs._grid->_rdimensions[dimension];

  if ( !rhs._grid->CheckerBoarded(dimension) ) {
    cbmask=0x3;
  }

  int ro  = rplane*rhs._grid->_ostride[dimension]; // base offset for start of plane 
  int lo  = lplane*lhs._grid->_ostride[dimension]; // base offset for start of plane 
  int o   = 0;                                     // relative offset to base within plane
  
#pragma omp parallel for collapse(2)
  for(int n=0;n<rhs._grid->_slice_nblock[dimension];n++){
    for(int b=0;b<rhs._grid->_slice_block[dimension];b++){
      
      int ocb=1<<lhs._grid->CheckerBoardFromOindex(o+b);
      if ( ocb&cbmask ) {
	permute(lhs._odata[lo+o+b],rhs._odata[ro+o+b],permute_type);
      }
      
    }
    o +=rhs._grid->_slice_stride[dimension];
  }
}

//////////////////////////////////////////////////////
// Local to node Cshift
//////////////////////////////////////////////////////
template<class vobj> void Cshift_local(Lattice<vobj>& ret,Lattice<vobj> &rhs,int dimension,int shift)
{
  int sshift[2];

  sshift[0] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,0);
  sshift[1] = rhs._grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,1);

  if ( sshift[0] == sshift[1] ) {
    Cshift_local(ret,rhs,dimension,shift,0x3);
  } else {
    Cshift_local(ret,rhs,dimension,shift,0x1);// if checkerboard is unfavourable take two passes
    Cshift_local(ret,rhs,dimension,shift,0x2);// both with block stride loop iteration
  }
}

template<class vobj> Lattice<vobj> Cshift_local(Lattice<vobj> &ret,Lattice<vobj> &rhs,int dimension,int shift,int cbmask)
{
  GridBase *grid = rhs._grid;
  int fd = grid->_fdimensions[dimension];
  int rd = grid->_rdimensions[dimension];
  int ld = grid->_ldimensions[dimension];
  int gd = grid->_gdimensions[dimension];

  // Map to always positive shift modulo global full dimension.
  shift = (shift+fd)%fd;

  ret.checkerboard = grid->CheckerBoardDestination(rhs.checkerboard,shift);
        
  // the permute type
  int permute_dim =grid->PermuteDim(dimension);
  int permute_type=grid->PermuteType(dimension);

  for(int x=0;x<rd;x++){       

    int o   = 0;
    int bo  = x * grid->_ostride[dimension];
    
    int cb= (cbmask==0x2)? 1 : 0;

    int sshift = grid->CheckerBoardShift(rhs.checkerboard,dimension,shift,cb);
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
