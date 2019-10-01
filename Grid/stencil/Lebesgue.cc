/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/stencil/Lebesgue.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/GridCore.h>
#include <algorithm>

NAMESPACE_BEGIN(Grid);

int LebesgueOrder::UseLebesgueOrder;
#ifdef KNL
std::vector<int> LebesgueOrder::Block({8,2,2,2});
#else
std::vector<int> LebesgueOrder::Block({2,2,2,2});
#endif
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

LebesgueOrder::LebesgueOrder(GridBase *_grid) 
{
  grid = _grid;
  if ( Block[0]==0) ZGraph();
  else if ( Block[1]==0) NoBlocking();
  else CartesianBlocking();

  if (0) {
    std::cout << "Thread Interleaving"<<std::endl;
    ThreadInterleave();
  } 
}
void LebesgueOrder::ThreadInterleave(void)
{
  Vector<IndexInteger> reorder = _LebesgueReorder;
  Vector<IndexInteger> throrder;
  int vol = _LebesgueReorder.size();
  int threads = GridThread::GetThreads();
  int blockbits=3;
  //  int blocklen = 8;
  //  int msk      = 0x7;

  for(int t=0;t<threads;t++){
    for(int ss=0;ss<vol;ss++){
      if ( ( ss >> blockbits) % threads == t ) { 
	throrder.push_back(reorder[ss]);
      }
    }
  }
  _LebesgueReorder = throrder;
}
void LebesgueOrder::NoBlocking(void) 
{
  std::cout<<GridLogDebug<<"Lexicographic : no cache blocking"<<std::endl;
  _LebesgueReorder.resize(0);
  for ( int s = 0 ; s!= grid->oSites();s++){
    _LebesgueReorder.push_back(s);
  }
}
void LebesgueOrder::CartesianBlocking(void) 
{
  _LebesgueReorder.resize(0);

  //    std::cout << GridLogDebug << " CartesianBlocking ";
  //    for(int d=0;d<Block.size();d++) std::cout <<Block[d]<<" ";
  //    std::cout<<std::endl; 

  IndexInteger ND = grid->_ndimension;

  assert(ND==4);
  assert(ND==Block.size());

  Coordinate dims(ND);
  Coordinate xo(ND,0);
  Coordinate xi(ND,0);

  for(IndexInteger mu=0;mu<ND;mu++){
    dims[mu] = grid->_rdimensions[mu];
  }

  IterateO(ND,ND-1,xo,xi,dims);
};

void LebesgueOrder::IterateO(int ND,int dim,
			     Coordinate & xo,
			     Coordinate & xi,
			     Coordinate &dims)
{
  for(xo[dim]=0;xo[dim]<dims[dim];xo[dim]+=Block[dim]){
    if ( dim > 0 ) {
      IterateO(ND,dim-1,xo,xi,dims);
    } else {
      IterateI(ND,ND-1,xo,xi,dims);
    }
  }
};

void LebesgueOrder::IterateI(int ND,
			     int dim,
			     Coordinate & xo,
			     Coordinate & xi,
			     Coordinate &dims)
{
  Coordinate x(ND);
  for(xi[dim]=0;xi[dim]<std::min(dims[dim]-xo[dim],Block[dim]);xi[dim]++){
    if ( dim > 0 ) {
      IterateI(ND,dim-1,xo,xi,dims);
    } else {
      for(int d=0;d<ND;d++){
	x[d]=xi[d]+xo[d];
	//	std::cout << x[d]<<" ";
      }
      //      std::cout << "\n";
      IndexInteger index;
      Lexicographic::IndexFromCoor(x,index,grid->_rdimensions);
      _LebesgueReorder.push_back(index);
    }
  }
}

void LebesgueOrder::ZGraph(void) 
{
  _LebesgueReorder.resize(0);

  std::cout << GridLogDebug << " Lebesgue order "<<std::endl;
  // Align up dimensions to power of two.
  const IndexInteger one=1;

  IndexInteger ND = grid->_ndimension;
  Coordinate dims(ND);
  Coordinate adims(ND);
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
  for(int bit=0;bit<32;bit++){
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

      assert(site < vol);
      _LebesgueReorder.push_back(site);
    }
  }
  assert( _LebesgueReorder.size() == vol );

  /*
    std::vector<int> coor(4);
    for(IndexInteger asite=0;asite<vol;asite++){
    grid->oCoorFromOindex (coor,_LebesgueReorder[asite]);
    std::cout << " site "<<asite << "->" << _LebesgueReorder[asite]<< " = ["
    << coor[0]<<","
    << coor[1]<<","
    << coor[2]<<","
    << coor[3]<<"]"
    <<std::endl;
    }
  */
}
NAMESPACE_END(Grid);

