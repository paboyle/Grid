/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonCompressor.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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
#ifndef  GRID_QCD_WILSON_COMPRESSOR_H
#define  GRID_QCD_WILSON_COMPRESSOR_H

NAMESPACE_BEGIN(Grid);

///////////////////////////////////////////////////////////////
// Wilson compressor will need FaceGather policies for:
// Periodic, Dirichlet, and partial Dirichlet for DWF
///////////////////////////////////////////////////////////////
const int dwf_compressor_depth=2;
#define DWF_COMPRESS
class FaceGatherPartialDWF
{
public:
#ifdef DWF_COMPRESS
  static int PartialCompressionFactor(GridBase *grid) {return grid->_fdimensions[0]/(2*dwf_compressor_depth);};
#else
  static int PartialCompressionFactor(GridBase *grid) { return 1;}
#endif
  template<class vobj,class cobj,class compressor>
  static void Gather_plane_simple (deviceVector<std::pair<int,int> >& table,
				   const Lattice<vobj> &rhs,
				   cobj *buffer,
				   compressor &compress,
				   int off,int so,int partial)
  {
    //DWF only hack: If a direction that is OFF node we use Partial Dirichlet
    //  Shrinks local and remote comms buffers
    GridBase *Grid = rhs.Grid();
    int Ls = Grid->_rdimensions[0];
#ifdef DWF_COMPRESS
    int depth=dwf_compressor_depth;
#else 
    int depth=Ls/2;
#endif
    std::pair<int,int> *table_v = & table[0];
    auto rhs_v = rhs.View(AcceleratorRead);
    int vol=table.size()/Ls;
    accelerator_forNB( idx,table.size(), vobj::Nsimd(), {
	Integer i=idx/Ls;
	Integer s=idx%Ls;
	Integer sc=depth+s-(Ls-depth);
	if(s<depth)     compress.Compress(buffer[off+i+s*vol],rhs_v[so+table_v[idx].second]);
	if(s>=Ls-depth) compress.Compress(buffer[off+i+sc*vol],rhs_v[so+table_v[idx].second]);
    });
    rhs_v.ViewClose();
  }
  template<class decompressor,class Decompression>
  static void DecompressFace(decompressor decompress,Decompression &dd)
  {
    auto Ls = dd.dims[0];
#ifdef DWF_COMPRESS
    int depth=dwf_compressor_depth;
#else
    int depth=Ls/2;
#endif    
    // Just pass in the Grid
    auto kp = dd.kernel_p;
    auto mp = dd.mpi_p;
    int size= dd.buffer_size;
    int vol= size/Ls;
    accelerator_forNB(o,size,1,{
	int idx=o/Ls;
	int   s=o%Ls;
	if ( s < depth ) {
	  int oo=s*vol+idx;
	  kp[o]=mp[oo];
	} else if ( s >= Ls-depth ) {
	  int sc = depth + s - (Ls-depth);
	  int oo=sc*vol+idx; 
	  kp[o]=mp[oo];
	} else {
	  kp[o] = Zero();//fill rest with zero if partial dirichlet
	}
    });
  }
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Need to gather *interior portions* for ALL s-slices in simd directions
  // Do the gather as need to treat SIMD lanes differently, and insert zeroes on receive side
  // Reorder the fifth dim to be s=Ls-1 , s=0, s=1,...,Ls-2.
  ////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj,class cobj,class compressor>
  static void Gather_plane_exchange(deviceVector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,
				    std::vector<cobj *> pointers,int dimension,int plane,int cbmask,
				    compressor &compress,int type,int partial)
  {
    GridBase *Grid = rhs.Grid();
    int Ls = Grid->_rdimensions[0];
#ifdef DWF_COMPRESS
    int depth=dwf_compressor_depth;
#else
    int depth = Ls/2;
#endif
    
    // insertion of zeroes...
    assert( (table.size()&0x1)==0);
    int num=table.size()/2;
    int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane
    
    auto rhs_v = rhs.View(AcceleratorRead);
    auto p0=&pointers[0][0];
    auto p1=&pointers[1][0];
    auto tp=&table[0];
    int nnum=num/Ls;
    accelerator_forNB(j, num, vobj::Nsimd(), {
	//  Reorders both local and remote comms buffers
	//  
	int s  = j % Ls;
	int sp1 = (s+depth)%Ls;  // peri incremented s slice
	
	int hxyz= j/Ls;

	int xyz0= hxyz*2; // xyzt part of coor
	int xyz1= hxyz*2+1;
	
	int jj= hxyz + sp1*nnum ; // 0,1,2,3 -> Ls-1 slice , 0-slice, 1-slice ....
	
	int kk0= xyz0*Ls + s ; // s=0 goes to s=1
	int kk1= xyz1*Ls + s ; // s=Ls-1 -> s=0
	compress.CompressExchange(p0[jj],p1[jj],
				  rhs_v[so+tp[kk0 ].second], // Same s, consecutive xyz sites
				  rhs_v[so+tp[kk1 ].second], 
				  type);
    });
    rhs_v.ViewClose();
  }
  // Merge routine is for SIMD faces
  template<class decompressor,class Merger>
  static void MergeFace(decompressor decompress,Merger &mm)
  {
    auto Ls = mm.dims[0];
#ifdef DWF_COMPRESS
    int depth=dwf_compressor_depth;
#else
    int depth = Ls/2;
#endif
    int  num= mm.buffer_size/2; // relate vol and Ls to buffer size
    auto mp = &mm.mpointer[0];
    auto vp0= &mm.vpointers[0][0]; // First arg is exchange first
    auto vp1= &mm.vpointers[1][0];
    auto type= mm.type;
    int nnum = num/Ls;
    accelerator_forNB(o,num,Merger::Nsimd,{

	int  s=o%Ls;
	int hxyz=o/Ls; // xyzt related component
	int xyz0=hxyz*2;
	int xyz1=hxyz*2+1;

	int sp = (s+depth)%Ls; 
	int jj= hxyz + sp*nnum ; // 0,1,2,3 -> Ls-1 slice , 0-slice, 1-slice ....

	int oo0= s+xyz0*Ls;
	int oo1= s+xyz1*Ls;

	// same ss0, ss1 pair goes to new layout
	decompress.Exchange(mp[oo0],mp[oo1],vp0[jj],vp1[jj],type);
      });
  }
};
class FaceGatherDWFMixedBCs
{
public:
#ifdef DWF_COMPRESS
  static int PartialCompressionFactor(GridBase *grid) {return grid->_fdimensions[0]/(2*dwf_compressor_depth);};
#else 
  static int PartialCompressionFactor(GridBase *grid) {return 1;}
#endif
  
  template<class vobj,class cobj,class compressor>
  static void Gather_plane_simple (deviceVector<std::pair<int,int> >& table,
					 const Lattice<vobj> &rhs,
					 cobj *buffer,
					 compressor &compress,
					 int off,int so,int partial)
  {
    //    std::cout << " face gather simple DWF partial "<<partial <<std::endl;
    if(partial) FaceGatherPartialDWF::Gather_plane_simple(table,rhs,buffer,compress,off,so,partial);
    else        FaceGatherSimple::Gather_plane_simple(table,rhs,buffer,compress,off,so,partial);
  }
  template<class vobj,class cobj,class compressor>
  static void Gather_plane_exchange(deviceVector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,
				    std::vector<cobj *> pointers,int dimension,int plane,int cbmask,
				    compressor &compress,int type,int partial)
  {
    //    std::cout << " face gather exch DWF partial "<<partial <<std::endl;
    if(partial) FaceGatherPartialDWF::Gather_plane_exchange(table,rhs,pointers,dimension, plane,cbmask,compress,type,partial);
    else        FaceGatherSimple::Gather_plane_exchange    (table,rhs,pointers,dimension, plane,cbmask,compress,type,partial);
  }
  template<class decompressor,class Merger>
  static void MergeFace(decompressor decompress,Merger &mm)
  {
    int partial = mm.partial;
    //    std::cout << " merge DWF partial "<<partial <<std::endl;
    if ( partial ) FaceGatherPartialDWF::MergeFace(decompress,mm);
    else           FaceGatherSimple::MergeFace(decompress,mm);
  }

  template<class decompressor,class Decompression>
  static void DecompressFace(decompressor decompress,Decompression &dd)
  {
    int partial = dd.partial;
    //    std::cout << " decompress DWF partial "<<partial <<std::endl;
    if ( partial ) FaceGatherPartialDWF::DecompressFace(decompress,dd);
    else           FaceGatherSimple::DecompressFace(decompress,dd);
  }
};

/////////////////////////////////////////////////////////////////////////////////////////////
// optimised versions supporting half precision too??? Deprecate
/////////////////////////////////////////////////////////////////////////////////////////////


//Could make FaceGather a template param, but then behaviour is runtime not compile time
template<class _HCspinor,class _Hspinor,class _Spinor, class projector>
class WilsonCompressorTemplate  : public FaceGatherDWFMixedBCs
//  : public FaceGatherSimple
{
public:
  
  int mu,dag;  

  void Point(int p) { mu=p; };

  WilsonCompressorTemplate(int _dag=0){
    dag = _dag;
  }

  typedef _Spinor         SiteSpinor;
  typedef _Hspinor     SiteHalfSpinor;
  typedef _HCspinor SiteHalfCommSpinor;
  typedef typename SiteHalfCommSpinor::vector_type vComplexLow;
  typedef typename SiteHalfSpinor::vector_type     vComplexHigh;
  constexpr static int Nw=sizeof(SiteHalfSpinor)/sizeof(vComplexHigh);

  accelerator_inline int CommDatumSize(void) const {
    return sizeof(SiteHalfCommSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  accelerator_inline void Compress(SiteHalfSpinor &buf,const SiteSpinor &in) const {
    typedef decltype(coalescedRead(buf)) sobj;
    sobj sp;
    auto sin = coalescedRead(in);
    projector::Proj(sp,sin,mu,dag);
    coalescedWrite(buf,sp);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  accelerator_inline void Exchange(SiteHalfSpinor &mp0,
				   SiteHalfSpinor &mp1,
				   const SiteHalfSpinor & vp0,
				   const SiteHalfSpinor & vp1,
				   Integer type) const {
#ifdef GRID_SIMT
    exchangeSIMT(mp0,mp1,vp0,vp1,type);
#else
    SiteHalfSpinor tmp1;
    SiteHalfSpinor tmp2;
    exchange(tmp1,tmp2,vp0,vp1,type);
    vstream(mp0,tmp1);
    vstream(mp1,tmp2);
#endif
  }
  

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  accelerator_inline void Decompress(SiteHalfSpinor &out,
				     SiteHalfSpinor &in) const {    
    out = in;
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  accelerator_inline void CompressExchange(SiteHalfSpinor &out0,
					   SiteHalfSpinor &out1,
					   const SiteSpinor &in0,
					   const SiteSpinor &in1,
					   Integer type) const
  {
#ifdef GRID_SIMT
    typedef SiteSpinor vobj;
    typedef SiteHalfSpinor hvobj;
    typedef decltype(coalescedRead(in0))    sobj;
    typedef decltype(coalescedRead(out0)) hsobj;

    constexpr unsigned int Nsimd = vobj::Nsimd();
    unsigned int mask = Nsimd >> (type + 1);
    int lane = acceleratorSIMTlane(Nsimd);
    int j0 = lane &(~mask); // inner coor zero
    int j1 = lane |(mask) ; // inner coor one
    const vobj *vp0 = &in0;
    const vobj *vp1 = &in1;
    const vobj *vp = (lane&mask) ? vp1:vp0;
    auto sa = coalescedRead(*vp,j0);
    auto sb = coalescedRead(*vp,j1);
    hsobj psa, psb;
    projector::Proj(psa,sa,mu,dag);
    projector::Proj(psb,sb,mu,dag);
    coalescedWrite(out0,psa);
    coalescedWrite(out1,psb);
#else
    SiteHalfSpinor temp1, temp2;
    SiteHalfSpinor temp3, temp4;
    projector::Proj(temp1,in0,mu,dag);
    projector::Proj(temp2,in1,mu,dag);
    exchange(temp3,temp4,temp1,temp2,type);
    vstream(out0,temp3);
    vstream(out1,temp4);
#endif
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  accelerator_inline bool DecompressionStep(void) const {
    return false;
  }

};

#define DECLARE_PROJ(Projector,Compressor,spProj)			\
  class Projector {							\
  public:								\
  template<class hsp,class fsp>						\
  static accelerator void Proj(hsp &result,const fsp &in,int mu,int dag){ \
    spProj(result,in);							\
  }									\
  };									\
  template<typename HCS,typename HS,typename S> using Compressor = WilsonCompressorTemplate<HCS,HS,S,Projector>;

DECLARE_PROJ(WilsonXpProjector,WilsonXpCompressor,spProjXp);
DECLARE_PROJ(WilsonYpProjector,WilsonYpCompressor,spProjYp);
DECLARE_PROJ(WilsonZpProjector,WilsonZpCompressor,spProjZp);
DECLARE_PROJ(WilsonTpProjector,WilsonTpCompressor,spProjTp);
DECLARE_PROJ(WilsonXmProjector,WilsonXmCompressor,spProjXm);
DECLARE_PROJ(WilsonYmProjector,WilsonYmCompressor,spProjYm);
DECLARE_PROJ(WilsonZmProjector,WilsonZmCompressor,spProjZm);
DECLARE_PROJ(WilsonTmProjector,WilsonTmCompressor,spProjTm);

class WilsonProjector {
public:
  template<class hsp,class fsp>
  static accelerator void Proj(hsp &result,const fsp &in,int mu,int dag){
    int mudag=dag? mu : (mu+Nd)%(2*Nd);
    switch(mudag) {
    case Xp:	spProjXp(result,in);	break;
    case Yp:	spProjYp(result,in);	break;
    case Zp:	spProjZp(result,in);	break;
    case Tp:	spProjTp(result,in);	break;
    case Xm:	spProjXm(result,in);	break;
    case Ym:	spProjYm(result,in);	break;
    case Zm:	spProjZm(result,in);	break;
    case Tm:	spProjTm(result,in);	break;
    default: 	assert(0);	        break;
    }
  }
};
template<typename HCS,typename HS,typename S> using WilsonCompressor = WilsonCompressorTemplate<HCS,HS,S,WilsonProjector>;

// Fast comms buffer manipulation which should inline right through (avoid direction
// dependent logic that prevents inlining
template<class vobj,class cobj,class Parameters>
class WilsonStencil : public CartesianStencil<vobj,cobj,Parameters> {
public:

  typedef CartesianStencil<vobj,cobj,Parameters> Base;
  typedef typename Base::View_type View_type;

  //  Vector<int> surface_list;
  WilsonStencil(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances,Parameters p)  
    : CartesianStencil<vobj,cobj,Parameters> (grid,npoints,checkerboard,directions,distances,p) 
  { 
    //    surface_list.resize(0);
    this->same_node.resize(npoints);
  };
  
  template < class compressor>
  void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
  {
    std::vector<std::vector<CommsRequest_t> > reqs;
    this->HaloExchangeOptGather(source,compress);
    // Asynchronous MPI calls multidirectional, Isend etc...
    // Non-overlapped directions within a thread. Asynchronous calls except MPI3, threaded up to comm threads ways.
    this->Communicate();
    this->CommsMerge(compress);
    this->CommsMergeSHM(compress);
  }
  
  template <class compressor>
  void HaloExchangeOptGather(const Lattice<vobj> &source,compressor &compress) 
  {
    this->Prepare();
    this->HaloGatherOpt(source,compress);
  }

  template <class compressor>
  void HaloGatherOpt(const Lattice<vobj> &source,compressor &compress)
  {
    // Strategy. Inherit types from Compressor.
    // Use types to select the write direction by directon compressor
    typedef typename compressor::SiteSpinor         SiteSpinor;
    typedef typename compressor::SiteHalfSpinor     SiteHalfSpinor;
    typedef typename compressor::SiteHalfCommSpinor SiteHalfCommSpinor;

    this->_grid->StencilBarrier();

    assert(source.Grid()==this->_grid);
    
    this->u_comm_offset=0;
      
    WilsonXpCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> XpCompress; 
    WilsonYpCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> YpCompress; 
    WilsonZpCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> ZpCompress; 
    WilsonTpCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> TpCompress;
    WilsonXmCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> XmCompress; 
    WilsonYmCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> YmCompress; 
    WilsonZmCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> ZmCompress; 
    WilsonTmCompressor<SiteHalfCommSpinor,SiteHalfSpinor,SiteSpinor> TmCompress;

    int dag = compress.dag;
    int face_idx=0;
#define vet_same_node(a,b) \
      { auto tmp = b;  }
    if ( dag ) { 
      vet_same_node(this->same_node[Xp],this->HaloGatherDir(source,XpCompress,Xp,face_idx));
      vet_same_node(this->same_node[Yp],this->HaloGatherDir(source,YpCompress,Yp,face_idx));
      vet_same_node(this->same_node[Zp],this->HaloGatherDir(source,ZpCompress,Zp,face_idx));
      vet_same_node(this->same_node[Tp],this->HaloGatherDir(source,TpCompress,Tp,face_idx));
      vet_same_node(this->same_node[Xm],this->HaloGatherDir(source,XmCompress,Xm,face_idx));
      vet_same_node(this->same_node[Ym],this->HaloGatherDir(source,YmCompress,Ym,face_idx));
      vet_same_node(this->same_node[Zm],this->HaloGatherDir(source,ZmCompress,Zm,face_idx));
      vet_same_node(this->same_node[Tm],this->HaloGatherDir(source,TmCompress,Tm,face_idx));
    } else {
      vet_same_node(this->same_node[Xp],this->HaloGatherDir(source,XmCompress,Xp,face_idx));
      vet_same_node(this->same_node[Yp],this->HaloGatherDir(source,YmCompress,Yp,face_idx));
      vet_same_node(this->same_node[Zp],this->HaloGatherDir(source,ZmCompress,Zp,face_idx));
      vet_same_node(this->same_node[Tp],this->HaloGatherDir(source,TmCompress,Tp,face_idx));
      vet_same_node(this->same_node[Xm],this->HaloGatherDir(source,XpCompress,Xm,face_idx));
      vet_same_node(this->same_node[Ym],this->HaloGatherDir(source,YpCompress,Ym,face_idx));
      vet_same_node(this->same_node[Zm],this->HaloGatherDir(source,ZpCompress,Zm,face_idx));
      vet_same_node(this->same_node[Tm],this->HaloGatherDir(source,TpCompress,Tm,face_idx));
    }
    this->face_table_computed=1;
    assert(this->u_comm_offset==this->_unified_buffer_size);
    accelerator_barrier();
  }

};

NAMESPACE_END(Grid);
#endif
