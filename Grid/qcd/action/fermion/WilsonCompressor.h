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

/////////////////////////////////////////////////////////////////////////////////////////////
// optimised versions supporting half precision too??? Deprecate
/////////////////////////////////////////////////////////////////////////////////////////////


//Could make FaceGather a template param, but then behaviour is runtime not compile time
template<class _HCspinor,class _Hspinor,class _Spinor, class projector>
class WilsonCompressorTemplate : public FaceGatherSimple
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
#ifdef NVLINK_GET
    this->_grid->StencilBarrier(); // He can now get mu local gather, I can get his
    // Synch shared memory on a single nodes; could use an asynchronous barrier here and defer check
    // Or issue barrier AFTER the DMA is running
#endif    
  }

};

NAMESPACE_END(Grid);
#endif
