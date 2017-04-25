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

namespace Grid {
namespace QCD {

/////////////////////////////////////////////////////////////////////////////////////////////
// optimised versions supporting half precision too
/////////////////////////////////////////////////////////////////////////////////////////////

template<class _HCspinor,class _Hspinor,class _Spinor, class projector,typename SFINAE = void >
class WilsonCompressorTemplate;


template<class _HCspinor,class _Hspinor,class _Spinor, class projector>
class WilsonCompressorTemplate< _HCspinor, _Hspinor, _Spinor, projector,
  typename std::enable_if<std::is_same<_HCspinor,_Hspinor>::value>::type >
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

  inline int CommDatumSize(void) {
    return sizeof(SiteHalfCommSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    projector::Proj(buf[o],in,mu,dag);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o) {    
    assert(0);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  inline void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(out0[j],out1[j],temp1,temp2,type);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return false; }

};

template<class _HCspinor,class _Hspinor,class _Spinor, class projector>
class WilsonCompressorTemplate< _HCspinor, _Hspinor, _Spinor, projector,
  typename std::enable_if<!std::is_same<_HCspinor,_Hspinor>::value>::type >
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

  inline int CommDatumSize(void) {
    return sizeof(SiteHalfCommSpinor);
  }

  /*****************************************************/
  /* Compress includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Compress(SiteHalfSpinor *buf,Integer o,const SiteSpinor &in) {
    SiteHalfSpinor hsp;
    SiteHalfCommSpinor *hbuf = (SiteHalfCommSpinor *)buf;
    projector::Proj(hsp,in,mu,dag);
    precisionChange((vComplexLow *)&hbuf[o],(vComplexHigh *)&hsp,Nw);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor *mp,
                       SiteHalfSpinor *vp0,
                       SiteHalfSpinor *vp1,
		       Integer type,Integer o){
    SiteHalfSpinor vt0,vt1;
    SiteHalfCommSpinor *vpp0 = (SiteHalfCommSpinor *)vp0;
    SiteHalfCommSpinor *vpp1 = (SiteHalfCommSpinor *)vp1;
    precisionChange((vComplexHigh *)&vt0,(vComplexLow *)&vpp0[o],Nw);
    precisionChange((vComplexHigh *)&vt1,(vComplexLow *)&vpp1[o],Nw);
    exchange(mp[2*o],mp[2*o+1],vt0,vt1,type);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor *out,
			 SiteHalfSpinor *in, Integer o){
    SiteHalfCommSpinor *hin=(SiteHalfCommSpinor *)in;
    precisionChange((vComplexHigh *)&out[o],(vComplexLow *)&hin[o],Nw);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  inline void CompressExchange(SiteHalfSpinor *out0,
			       SiteHalfSpinor *out1,
			       const SiteSpinor *in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    SiteHalfCommSpinor *hout0 = (SiteHalfCommSpinor *)out0;
    SiteHalfCommSpinor *hout1 = (SiteHalfCommSpinor *)out1;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(temp3,temp4,temp1,temp2,type);
    precisionChange((vComplexLow *)&hout0[j],(vComplexHigh *)&temp3,Nw);
    precisionChange((vComplexLow *)&hout1[j],(vComplexHigh *)&temp4,Nw);
  }

  /*****************************************************/
  /* Pass the info to the stencil */
  /*****************************************************/
  inline bool DecompressionStep(void) { return true; }

};

#define DECLARE_PROJ(Projector,Compressor,spProj)			\
  class Projector {							\
  public:								\
    template<class hsp,class fsp>					\
    static void Proj(hsp &result,const fsp &in,int mu,int dag){			\
      spProj(result,in);						\
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
  static void Proj(hsp &result,const fsp &in,int mu,int dag){
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
template<class vobj,class cobj>
class WilsonStencil : public CartesianStencil<vobj,cobj> {
public:

  typedef CartesianCommunicator::CommsRequest_t CommsRequest_t;

  std::vector<int> same_node;
  std::vector<int> surface_list;

  WilsonStencil(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances)  
    : CartesianStencil<vobj,cobj> (grid,npoints,checkerboard,directions,distances) ,
    same_node(npoints)
  { 
    surface_list.resize(0);
  };

  void BuildSurfaceList(int Ls,int vol4){

    // find same node for SHM
    // Here we know the distance is 1 for WilsonStencil
    for(int point=0;point<this->_npoints;point++){
      same_node[point] = this->SameNode(point);
      //      std::cout << " dir " <<point<<" same_node " <<same_node[point]<<std::endl;
    }
    
    for(int site = 0 ;site< vol4;site++){
      int local = 1;
      for(int point=0;point<this->_npoints;point++){
	if( (!this->GetNodeLocal(site*Ls,point)) && (!same_node[point]) ){ 
	  local = 0;
	}
      }
      if(local == 0) { 
	surface_list.push_back(site);
      }
    }
  }

  template < class compressor>
  void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
  {
    std::vector<std::vector<CommsRequest_t> > reqs;
    this->HaloExchangeOptGather(source,compress);
    this->CommunicateBegin(reqs);
    this->CommunicateComplete(reqs);
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

    assert(source._grid==this->_grid);
    this->halogtime-=usecond();
    
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
    if ( dag ) { 
      //	std::cout << " Optimised Dagger compress " <<std::endl;
      assert(same_node[Xp]==this->HaloGatherDir(source,XpCompress,Xp,face_idx));
      assert(same_node[Yp]==this->HaloGatherDir(source,YpCompress,Yp,face_idx));
      assert(same_node[Zp]==this->HaloGatherDir(source,ZpCompress,Zp,face_idx));
      assert(same_node[Tp]==this->HaloGatherDir(source,TpCompress,Tp,face_idx));
      assert(same_node[Xm]==this->HaloGatherDir(source,XmCompress,Xm,face_idx));
      assert(same_node[Ym]==this->HaloGatherDir(source,YmCompress,Ym,face_idx));
      assert(same_node[Zm]==this->HaloGatherDir(source,ZmCompress,Zm,face_idx));
      assert(same_node[Tm]==this->HaloGatherDir(source,TmCompress,Tm,face_idx));
    } else {
      assert(same_node[Xp]==this->HaloGatherDir(source,XmCompress,Xp,face_idx));
      assert(same_node[Yp]==this->HaloGatherDir(source,YmCompress,Yp,face_idx));
      assert(same_node[Zp]==this->HaloGatherDir(source,ZmCompress,Zp,face_idx));
      assert(same_node[Tp]==this->HaloGatherDir(source,TmCompress,Tp,face_idx));
      assert(same_node[Xm]==this->HaloGatherDir(source,XpCompress,Xm,face_idx));
      assert(same_node[Ym]==this->HaloGatherDir(source,YpCompress,Ym,face_idx));
      assert(same_node[Zm]==this->HaloGatherDir(source,ZpCompress,Zm,face_idx));
      assert(same_node[Tm]==this->HaloGatherDir(source,TpCompress,Tm,face_idx));
    }
    this->face_table_computed=1;
    assert(this->u_comm_offset==this->_unified_buffer_size);
    this->halogtime+=usecond();
  }

 };

}} // namespace close
#endif
