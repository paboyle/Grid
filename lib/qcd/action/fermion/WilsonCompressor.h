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
  inline void Compress(SiteHalfSpinor * __restrict__ buf,Integer o,const SiteSpinor &in) {
    SiteHalfSpinor tmp;
    projector::Proj(tmp,in,mu,dag);
    vstream(buf[o],tmp);
  }

  /*****************************************************/
  /* Exchange includes precision change if mpi data is not same */
  /*****************************************************/
  inline void Exchange(SiteHalfSpinor * __restrict__ mp,
                       const SiteHalfSpinor * __restrict__ vp0,
                       const SiteHalfSpinor * __restrict__ vp1,
		       Integer type,Integer o){
    SiteHalfSpinor tmp1;
    SiteHalfSpinor tmp2;
    exchange(tmp1,tmp2,vp0[o],vp1[o],type);
    vstream(mp[2*o  ],tmp1);
    vstream(mp[2*o+1],tmp2);
  }

  /*****************************************************/
  /* Have a decompression step if mpi data is not same */
  /*****************************************************/
  inline void Decompress(SiteHalfSpinor * __restrict__ out,
			 SiteHalfSpinor * __restrict__ in, Integer o) {    
    assert(0);
  }

  /*****************************************************/
  /* Compress Exchange                                 */
  /*****************************************************/
  inline void CompressExchange(SiteHalfSpinor * __restrict__ out0,
			       SiteHalfSpinor * __restrict__ out1,
			       const SiteSpinor * __restrict__ in,
			       Integer j,Integer k, Integer m,Integer type){
    SiteHalfSpinor temp1, temp2,temp3,temp4;
    projector::Proj(temp1,in[k],mu,dag);
    projector::Proj(temp2,in[m],mu,dag);
    exchange(temp3,temp4,temp1,temp2,type);
    vstream(out0[j],temp3);
    vstream(out1[j],temp4);
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
  double timer0;
  double timer1;
  double timer2;
  double timer3;
  double timer4;
  double timer5;
  double timer6;
  uint64_t callsi;
  void ZeroCountersi(void)
  {
    timer0=0;
    timer1=0;
    timer2=0;
    timer3=0;
    timer4=0;
    timer5=0;
    timer6=0;
    callsi=0;
  }
  void Reporti(int calls)
  {
    if ( timer0 ) std::cout << GridLogMessage << " timer0 (HaloGatherOpt) " <<timer0/calls <<std::endl;
    if ( timer1 ) std::cout << GridLogMessage << " timer1 (Communicate)   " <<timer1/calls <<std::endl;
    if ( timer2 ) std::cout << GridLogMessage << " timer2 (CommsMerge )   " <<timer2/calls <<std::endl;
    if ( timer3 ) std::cout << GridLogMessage << " timer3 (commsMergeShm) " <<timer3/calls <<std::endl;
    if ( timer4 ) std::cout << GridLogMessage << " timer4 " <<timer4 <<std::endl;
  }

  WilsonStencil(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances)  
    : CartesianStencil<vobj,cobj> (grid,npoints,checkerboard,directions,distances) 
  { 
    ZeroCountersi();
  };


  template < class compressor>
  void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
  {
    std::vector<std::vector<CommsRequest_t> > reqs;
    this->HaloExchangeOptGather(source,compress);
    double t1=usecond();
    // Asynchronous MPI calls multidirectional, Isend etc...
    //    this->CommunicateBegin(reqs);
    //    this->CommunicateComplete(reqs);
    // Non-overlapped directions within a thread. Asynchronous calls except MPI3, threaded up to comm threads ways.
    this->Communicate();
    double t2=usecond(); timer1 += t2-t1;
    this->CommsMerge(compress);
    double t3=usecond(); timer2 += t3-t2;
    this->CommsMergeSHM(compress);
    double t4=usecond(); timer3 += t4-t3;
  }
  
  template <class compressor>
  void HaloExchangeOptGather(const Lattice<vobj> &source,compressor &compress) 
  {
    this->Prepare();
    double t0=usecond();
    this->HaloGatherOpt(source,compress);
    double t1=usecond();
    timer0 += t1-t0;
    callsi++;
  }

  template <class compressor>
  void HaloGatherOpt(const Lattice<vobj> &source,compressor &compress)
  {
    // Strategy. Inherit types from Compressor.
    // Use types to select the write direction by directon compressor
    typedef typename compressor::SiteSpinor         SiteSpinor;
    typedef typename compressor::SiteHalfSpinor     SiteHalfSpinor;
    typedef typename compressor::SiteHalfCommSpinor SiteHalfCommSpinor;

    this->mpi3synctime_g-=usecond();
    this->_grid->StencilBarrier();
    this->mpi3synctime_g+=usecond();

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
      assert(this->same_node[Xp]==this->HaloGatherDir(source,XpCompress,Xp,face_idx));
      assert(this->same_node[Yp]==this->HaloGatherDir(source,YpCompress,Yp,face_idx));
      assert(this->same_node[Zp]==this->HaloGatherDir(source,ZpCompress,Zp,face_idx));
      assert(this->same_node[Tp]==this->HaloGatherDir(source,TpCompress,Tp,face_idx));
      assert(this->same_node[Xm]==this->HaloGatherDir(source,XmCompress,Xm,face_idx));
      assert(this->same_node[Ym]==this->HaloGatherDir(source,YmCompress,Ym,face_idx));
      assert(this->same_node[Zm]==this->HaloGatherDir(source,ZmCompress,Zm,face_idx));
      assert(this->same_node[Tm]==this->HaloGatherDir(source,TmCompress,Tm,face_idx));
    } else {
      assert(this->same_node[Xp]==this->HaloGatherDir(source,XmCompress,Xp,face_idx));
      assert(this->same_node[Yp]==this->HaloGatherDir(source,YmCompress,Yp,face_idx));
      assert(this->same_node[Zp]==this->HaloGatherDir(source,ZmCompress,Zp,face_idx));
      assert(this->same_node[Tp]==this->HaloGatherDir(source,TmCompress,Tp,face_idx));
      assert(this->same_node[Xm]==this->HaloGatherDir(source,XpCompress,Xm,face_idx));
      assert(this->same_node[Ym]==this->HaloGatherDir(source,YpCompress,Ym,face_idx));
      assert(this->same_node[Zm]==this->HaloGatherDir(source,ZpCompress,Zm,face_idx));
      assert(this->same_node[Tm]==this->HaloGatherDir(source,TpCompress,Tm,face_idx));
    }
    this->face_table_computed=1;
    assert(this->u_comm_offset==this->_unified_buffer_size);
    this->halogtime+=usecond();
  }

 };

}} // namespace close
#endif
