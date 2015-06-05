#include <Grid.h>

namespace Grid {

int LebesgueOrder::UseLebesgueOrder;

LebesgueOrder::IndexInteger LebesgueOrder::alignup(IndexInteger n){
  n--;           // 1000 0011 --> 1000 0010
  n |= n >> 1;   // 1000 0010 | 0100 0001 = 1100 0011
  n |= n >> 2;   // 1100 0011 | 0011 0000 = 1111 0011
  n |= n >> 4;   // 1111 0011 | 0000 1111 = 1111 1111
  n |= n >> 8;   // ... (At this point all bits are 1, so further bitwise-or
  n |= n >> 16;  //      operations produce no effect.)
  n++;           // 1111 1111 --> 1 0000 0000
  return n;
};

LebesgueOrder::LebesgueOrder(GridBase *grid) 
{
  _LebesgueReorder.resize(0);
  
  // Align up dimensions to power of two.
  const IndexInteger one=1;
  IndexInteger ND = grid->_ndimension;
  std::vector<IndexInteger> dims(ND);
  std::vector<IndexInteger> adims(ND);
  std::vector<std::vector<IndexInteger> > bitlist(ND);
  
  for(IndexInteger mu=0;mu<ND;mu++){
    dims[mu] = grid->_rdimensions[mu];
    assert ( dims[mu] != 0 );
    adims[mu] = alignup(dims[mu]);
  }
  
  // List which bits of padded volume coordinate contribute; this strategy 
  // i) avoids recursion 
  // ii) has loop lengths at most the width of a 32 bit word.
  int sitebit=0;
  int split=2;
  for(int mu=0;mu<ND;mu++){   // mu 0 takes bit 0; mu 1 takes bit 1 etc...
    for(int bit=0;bit<split;bit++){
      IndexInteger mask = one<<bit;
      if ( mask&(adims[mu]-1) ){
	bitlist[mu].push_back(sitebit);
	sitebit++;
      }
    }
  }
  for(int bit=split;bit<32;bit++){
    IndexInteger mask = one<<bit;
    for(int mu=0;mu<ND;mu++){   // mu 0 takes bit 0; mu 1 takes bit 1 etc...
      if ( mask&(adims[mu]-1) ){
	bitlist[mu].push_back(sitebit);
	sitebit++;
      }
    }
  }
  
  // Work out padded and unpadded volumes
  IndexInteger avol = 1;
  for(int mu=0;mu<ND;mu++) avol = avol * adims[mu];
  
  IndexInteger vol = 1;
  for(int mu=0;mu<ND;mu++) vol = vol * dims[mu];
  
  // Loop over padded volume, following Lebesgue curve
  // We interleave the bits from sequential "mu".
  std::vector<IndexInteger> ax(ND);
  
  for(IndexInteger asite=0;asite<avol;asite++){
    
    // Start with zero and collect bits
    for(int mu=0;mu<ND;mu++) ax[mu] = 0;
    
    int contained = 1;
    for(int mu=0;mu<ND;mu++){
      
      // Build the coordinate on the aligned volume
      for(int bit=0;bit<bitlist[mu].size();bit++){
	int sbit=bitlist[mu][bit];
	
	if(asite&(one<<sbit)){
	  ax[mu]|=one<<bit;
	}
      }
      
      // Is it contained in original box
      if ( ax[mu]>dims[mu]-1 ) contained = 0;
      
    }
    
    if ( contained ) {
      int site = ax[0]
	+        dims[0]*ax[1]
	+dims[0]*dims[1]*ax[2]
	+dims[0]*dims[1]*dims[2]*ax[3];
      
      _LebesgueReorder.push_back(site);
    }
	}
  assert( _LebesgueReorder.size() == vol );
}
}
